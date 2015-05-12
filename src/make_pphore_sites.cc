//
// file make_pphore_sites
// David Cosgrove
// AstraZeneca
// 30th March 2006
//
// This file contains the functions to build a vector of vector of
// SinglePPhoreSites, 1 vector per conformation in the given OEMol.

#include <iostream>
#include <oechem.h>

#include "gtpl_defs.H"
#include "points_on_sphere.H"
#include "stddefs.H"
#include "OverlayTrans.H"
#include "PharmPoint.H"
#include "SinglePPhoreSite.H"
#include "SMARTSExceptions.H"

#include <boost/lexical_cast.hpp>

using namespace std;
using namespace OEChem;
using namespace OESystem;
using namespace OEPlatform;

// take the _Conf<nn> off the end of the molecule name - in oemol_conf_name_fns.cc
string root_mol_name( const string &full_mol_name );

namespace DACLIB {

// for find_hphobes_itmoc, it's the minimum hydrophobicity value for a group to
// be counted as a hydrophobe 
static const float H_MIN = 8.7;

void rotate_about_axis( const float axis[3] , float angle , float vec[3] );
void rotate_about_axis( const double axis[3] , double angle , double vec[3] );
double angle_about_axis( const float axis[3] , const float vec1[3] ,
                         const float vec2[3] );
double angle_about_axis( const double axis[3] , const double vec1[3] ,
                         const double vec2[3] );

// ****************************************************************************
// return the type number of the ITMOC points
int get_itmoc_type_num( PharmPoint &pharm_points ) {

  int i = 0;
  map<string,vector<string> >::iterator q;
  for( i = 0 , q = pharm_points.points_defs().begin() ;
       q != pharm_points.points_defs().end() ; ++i , ++q ) {
    if( pharm_points.hphobes_itmoc() &&
        pharm_points.itmoc_label() == q->first ) {
      break;
    }
    if( pharm_points.hphobes_itmoc_alo() &&
        pharm_points.itmoc_alo_label() == q->first ) {
      break;
    }
  }

  return i;

}

// ****************************************************************************
// throws DACLIB::SMARTSDefnError exceptions if there's any errors with the SMARTS.
void build_smarts_object( const string &smarts_label ,
                          const vector<pair<string,string> > &input_smarts ,
                          vector<pair<string,string> > &smarts_sub_defn ,
                          map<string , OESubSearch *> &smts_objs ) {

  if( input_smarts.empty() ) {
    // bail soonest
    throw( DACLIB::SMARTSDefnError( "No SMARTS definitions." ) );
  }

  // find the SMARTS string corresponding to the label
  vector<pair<string,string> >::const_iterator p;
  string short_smarts;
  for( p = input_smarts.begin() ; p != input_smarts.end() ; ++p ) {
    if( p->first == smarts_label ) {
      short_smarts = p->second;
      break;
    }
  }

  // do any vector binding expansion
  string exp_smarts = short_smarts;
  if( string::npos != exp_smarts.find( "$" ) &&
      !OESmartsLexReplace( exp_smarts , smarts_sub_defn ) ) {
    // We shouldn't have got this far without all the necessary sub-definitions
    // so just bail
    ostringstream oss;
    oss << "AWOOGA - something screwed with the SMARTS sub-definitions for "
        << endl << smarts_label << " : " << short_smarts << endl
        << "in DACLIB::build_smarts_object()" << endl;
    throw( DACLIB::SMARTSDefnError( oss.str().c_str() ) );
  }

  // don't allow re-ordering of SMARTS for efficiency - we might be relying on
  // the hit atoms coming out in the SMARTS order.
  OESubSearch *smarts_obj = new OESubSearch( exp_smarts.c_str() , false );
  if( !*smarts_obj ) {
    ostringstream oss;
    oss << "Error parsing SMARTS " << endl
        << smarts_label << " : " << short_smarts
        << " in DACLIB::build_smarts_objects()." << endl
        << "It expanded to " << exp_smarts << endl;
    throw( DACLIB::SMARTSDefnError( oss.str().c_str() ) );
  }

  smts_objs.insert( make_pair( smarts_label , smarts_obj ) );

}

// *********************************************************************
bool delocalised_o_or_n( OEAtomBase *target_atom ,
                         OEAtomBase *&twiddle_base ) {

  if( target_atom->GetAtomicNum() != OEElemNo::N &&
      target_atom->GetAtomicNum() != OEElemNo::O ) {
    return false; // obviously, it's not
  }

  // check for multiple bonds
  OEIter<OEBondBase> bond;
  for( bond = target_atom->GetBonds() ; bond ; ++bond ) {
    if( bond->GetOrder() > 1 || bond->IsAromatic() ) {
      if( bond->GetBgnIdx() == target_atom->GetIdx() )
        twiddle_base = bond->GetEnd();
      else
        twiddle_base = bond->GetBgn();
      return true;
    }
  }

  // check for N or O attached to unsaturated C
  twiddle_base = 0;
  OEIter<OEAtomBase> last_base = target_atom->GetAtoms( OEIsHeavy() );
  for( ; last_base ; ++last_base ) {
    if( ( OEElemNo::N == target_atom->GetAtomicNum() ||
          OEElemNo::O == target_atom->GetAtomicNum() ) ) {
      if( OEElemNo::C == last_base->GetAtomicNum() ) {
        twiddle_base = last_base;
        if( last_base->IsAromatic() ) {
          return true;
	}
        for( bond = last_base->GetBonds() ; bond ; ++bond ) {
          if( bond->GetOrder() > 1 ) {
            return true;
	  }
	}
      } else {
        // anything else is not delocalised, by my reckoning
        OEIter<OEAtomBase> tb = target_atom->GetAtoms( OEIsHeavy() );
        twiddle_base = tb;
        return false;
      }
    }
  }

  return false;

}

// ************************************************************************
// determine whether the given atom in the given molecule is twiddlable - it's
// not if it's connected to more than 1 heavy atom or by a multiple bond or if
// it's an N attached to an unsaturated C, such as amide and aniline or O likewise.
bool is_atom_twiddlable( OEMol &mol , unsigned int at_num ,
                         OEAtomBase *&twiddle_base ) {

  // the atom number passed in is the result of GetIdx(), not a sequence number
  // so need to find the atom the hard way
  OEIter<OEAtomBase> atom;
  atom = mol.GetAtoms( OEHasAtomIdx( at_num ) );
  if( !atom || !atom->GetHvyDegree() ) {
    return false; // didn't find the atom, or it has no n'bours
  }
  OEAtomBase *target_atom = atom;

  // any double bonds off the atom?
  OEIter<OEBondBase> bond = target_atom->GetBonds();
  for( ; bond ; ++bond ) {
    if( bond->GetOrder() != 1 ) {
      return false;
    }
  }

  if( target_atom->GetHvyDegree() > 1 ) {
    return false; // more than 1 heavy atom
  }

#ifdef NOTYET
  cout << "target atom : " << target_atom->GetIdx() << " : " << target_atom->GetAtomicNum() << endl;
#endif
  if( delocalised_o_or_n( target_atom , twiddle_base ) ) {
    return false;
  } else {
    // if it isn't delocalised o or n, and twiddle_base has been set non-zero,
    // we're good to go
    if( twiddle_base ) {
      return true;
    }
  }

  // final catch-all
  atom = target_atom->GetAtoms();
  twiddle_base = atom;

  return true;

}

// ************************************************************************
// determine whether the given atom in the given molecule is flippable. This
// is so if it's a phenolic O/S or an enamine.
bool is_atom_flippable( OEMol &mol , unsigned int at_num ,
                        OEAtomBase *&twiddle_base ) {

  OEIter<OEAtomBase> atom;
  atom = mol.GetAtoms( OEHasAtomIdx( at_num ) );
  if( !atom || !atom->GetHvyDegree() ) {
    return false; // didn't find the atom, or is has no h'bours
  }
  OEAtomBase *target_atom = atom;

  if( OEElemNo::N == target_atom->GetAtomicNum() &&
      1 == target_atom->GetTotalHCount() &&
      1 == target_atom->GetHvyDegree() ) {
    OEIter<OEAtomBase> tb = target_atom->GetAtoms( OEIsHeavy() );
    twiddle_base = tb;
    return true;
  }

  if( ( OEElemNo::O == target_atom->GetAtomicNum() ||
        OEElemNo::S == target_atom->GetAtomicNum() ) &&
      1 == target_atom->GetTotalHCount() && 1 == target_atom->GetHvyDegree() ) {
    OEIter<OEAtomBase> tb = target_atom->GetAtoms( OEIsHeavy() );
    twiddle_base = tb;
    return twiddle_base->IsAromatic();
  }

  return false;

}

// ************************************************************************
// make a SinglePPhoreSite at the centre of the given coords and put it in the list.
void create_pphore_site( double *at_cds , const int *site_atoms ,
                         int num_site_atoms , int type_code ,
                         const string &type_string , const string &label ,
                         const string &mol_name ,
                         vector<SinglePPhoreSite *> &pharm_sites ) {

  int      i;

  double    cds[3] = { 0.0 , 0.0 , 0.0 };
  double    dir[3] = { 0.0 , 0.0 , 0.0 };

  double *pat_cds = at_cds;
  for( i = 0 ; i < num_site_atoms ; i++ ) {
    cds[0] += *(pat_cds++);
    cds[1] += *(pat_cds++);
    cds[2] += *(pat_cds++);
  }

  cds[0] /= float( num_site_atoms );
  cds[1] /= float( num_site_atoms );
  cds[2] /= float( num_site_atoms );

  try {

    SinglePPhoreSite *ps =
        new SinglePPhoreSite( cds , dir , type_code , type_string , label ,
                              false , num_site_atoms , site_atoms , mol_name );
    pharm_sites.push_back( ps );
    ps->set_selected( false );

  } catch( int max_atoms ) {
    // need to split it up, really, but just print an error for now and the point
    // won't be made.
    cerr << "site " << label << " in molecule " << mol_name
         << " has too many atoms : got " << num_site_atoms
         << " but can't have more than " << MAX_SITE_ATOMS << endl;
  }

}

// ************************************************************************
// create a SinglePPhoreSite with a normal given the picked atoms
void create_pphore_site_with_normal( double *at_cds , const int *site_atoms ,
                                     int num_site_atoms , int type_code ,
                                     const string &type_string ,
                                     const string &label ,
                                     const string &mol_name ,
                                     vector<SinglePPhoreSite *> &pharm_sites ) {

  int      i;

  create_pphore_site( at_cds , site_atoms , num_site_atoms , type_code ,
                      type_string , label , root_mol_name( mol_name ) ,
                      pharm_sites );

  SinglePPhoreSite *ps = pharm_sites.back();

  const double *c_cds = ps->coords();
  double vec1[3] , vec2[3] , this_norm[3] , norm[3] = { 0.0 , 0.0 , 0.0 };
  double dat_cds1[3] , dat_cds2[3];

  double *pat_cds = at_cds;

  for( i = 0 ; i < num_site_atoms - 1 ; i++ ) {

    dat_cds1[0] = *(pat_cds++);
    dat_cds1[1] = *(pat_cds++);
    dat_cds1[2] = *(pat_cds++);

    dat_cds2[0] = pat_cds[0];
    dat_cds2[1] = pat_cds[1];
    dat_cds2[2] = pat_cds[2];

    DACLIB::join_vector( c_cds , dat_cds1 , vec1 );
    DACLIB::join_vector( c_cds , dat_cds2 , vec2 );
    DACLIB::norm_cross_product( vec1 , vec2 , this_norm );
    norm[0] += this_norm[0];
    norm[1] += this_norm[1];
    norm[2] += this_norm[2];

  }
  
  dat_cds1[0] = pat_cds[0];
  dat_cds1[1] = pat_cds[1];
  dat_cds1[2] = pat_cds[2];

  dat_cds2[0] = at_cds[0];
  dat_cds2[1] = at_cds[1];
  dat_cds2[2] = at_cds[2];

  DACLIB::join_vector( c_cds , dat_cds1 , vec1 );
  DACLIB::join_vector( c_cds , dat_cds2 , vec2 );
  DACLIB::norm_cross_product( vec1 , vec2 , this_norm );
  
  norm[0] += this_norm[0];
  norm[1] += this_norm[1];
  norm[2] += this_norm[2];

  norm[0] /= double( num_site_atoms );
  norm[1] /= double( num_site_atoms );
  norm[2] /= double( num_site_atoms );

  double l = DACLIB::length( norm );
  norm[0] /= l; norm[1] /= l; norm[2] /= l;

  ps->set_direction( norm , GtplDefs::RING_NORMAL );

} 

// ***************************************************************************
// build a tetrahedron such that one of its legs points along the target
// vector. Requires that otd can hold 15 floats - the ends of the 4 legs and
// the centre.
void make_oriented_tetrahedron( float target_vec[3] , float *otd ) {

  // centred on origin, legs pointing out from there.
  static const float td[15] = { 0.0F , 0.0F , 0.0F ,
                                0.0 , 0.9246 , -0.3335 ,
                                0.8163 , -0.4713 , -0.3335 ,
                                -0.8163 , -0.4713 , -0.3335 ,
                                0.0 , 0.0 , 1.0 };

  copy( td , td + 15 , otd );

  // if it's already along the z axis, nowt to do
  if( fabs( target_vec[0] ) < 1e-10 && fabs( target_vec[1] ) < 1e-10 ) {
    if( fabs( target_vec[2] - 1.0 ) < 1e-10 )
      return;
    else if( fabs( target_vec[2] + 1.0 ) < 1e-10 ) {
      // pointing along -z, so reflect
      otd[2] *= -1.0;
      otd[5] *= -1.0;
      otd[8] *= -1.0;
      otd[11] *= -1.0;
      otd[14] *= -1.0;
      return;
    }
  }

  // find the cross-product between the z axis and target_vec and the angle
  // between them, which is the angle that the z axis needs to be rotated by
  // about the cross-product to get to where we want it. Doing it this way
  // involves 2 cross-products, but it gets the sign of the angle correct.
  float cp[3];
  static const float *z_axis = td + 12;

  DACLIB::cross_product( target_vec , z_axis , cp );

  float a = -DACLIB::angle( target_vec , 1.0F , z_axis , 1.0F );

  DACLIB::rotate_about_axis( cp , a , otd + 3 );
  DACLIB::rotate_about_axis( cp , a , otd + 6 );
  DACLIB::rotate_about_axis( cp , a , otd + 9 );
  DACLIB::rotate_about_axis( cp , a , otd + 12 );

}

// ************************************************************************
bool get_attachment_vector( OEAtomBase *target_atom , OEMolBase &mol ,
                            float att_vec[3] ) {

  OEIter<OEAtomBase> atom = target_atom->GetAtoms( OEIsHeavy() );
  OEAtomBase *x = atom;
  if( !x )
    return false;

  float at_cds[3] , x_cds[3];
  mol.GetCoords( target_atom , at_cds );
  mol.GetCoords( x , x_cds );

  DACLIB::join_vector( x_cds , at_cds , att_vec );
  DACLIB::normalise( att_vec );

  return true;

}

// ************************************************************************
// count the number of explicit H atoms that are not at the same coords
// as the target atom (which is where OEChem puts them, sometimes)
int count_exp_h_atoms_at_dist( OEMolBase &mol , OEAtomBase *target_atom ) {

  int num_h = 0;
  OEIter<OEAtomBase> h = target_atom->GetAtoms( OEHasAtomicNum( OEElemNo::H ) );
  for( ; h ; ++h ) {
    if( OEGetDistance( mol , target_atom , h ) > 0.0 ) {
      ++num_h;
    }
  }

  return num_h;

}

// ************************************************************************
// add the required number of h vectors to the site, based on the direction
// of the target->explicit H atom bonds. We're assuming there'll be enough,
// as checked by count_exp_h_atoms_at_dist.
// Not using this for now, as Omega makes pyramidal aniline hydrogens, which
// are too wrong to be tolerated. But leaving it in in case this is ever
// fixed.
void add_h_vector_using_exp_h( OEMolBase &mol , OEAtomBase *target_atom ,
                               int num_h_to_add , SinglePPhoreSite &site  ) {

  float t_cds[3] , v[3];
  mol.GetCoords( target_atom , t_cds );
  OEIter<OEAtomBase> h = target_atom->GetAtoms( OEHasAtomicNum( OEElemNo::H ) );
  int num_h_added = 0;
  for( ; h && num_h_added < num_h_to_add ; ++h , ++num_h_added ) {
    float h_cds[3];
    mol.GetCoords( h , h_cds );
    DACLIB::join_vector( t_cds , h_cds , v );
    DACLIB::normalise( v );
    site.set_direction( v , GtplDefs::H_VECTOR );
  }

}

// ************************************************************************
void add_tet_h_vector_to_o_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                 SinglePPhoreSite &site ) {

  float att_vec[3];
  if( !get_attachment_vector( target_atom , mol , att_vec ) )
    return;

  float td[15];
  make_oriented_tetrahedron( att_vec , td );

  float d[3];
  d[0] = -td[3];
  d[1] = -td[4];
  d[2] = -td[5];
  site.set_direction( d , GtplDefs::H_VECTOR );

}

// ************************************************************************
void add_planar_h_vector_to_o_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                    OEAtomBase *base_atom ,
                                    SinglePPhoreSite &site ) {

  // first, add a direction
  add_tet_h_vector_to_o_site( mol , target_atom , site );

  // need to line the direction up with the plane of the base_atom bonds
  const double *dir = 0;
  int dir_num = -1;
  for( int i = 0 , is = site.get_num_dirs() ; i < is ; ++i ) {
    if( GtplDefs::H_VECTOR == site.direction_type( i ) ) {
      dir = site.direction( i );
      dir_num = i;
      break;
    }
  }
  if( !dir ) {
    return; // something went wrong.
  }

  OEIter<OEAtomBase> atom = base_atom->GetAtoms( OEIsHeavy() );
  if( atom->GetIdx() == target_atom->GetIdx() || 1 == atom->GetAtomicNum() ) {
    ++atom;
  }
  if( !atom ) {
    return; // wierd, though, what?
  }
  OEAtomBase *x_atom = atom;

  float b_cds[3] , t_cds[3] , x_cds[3];
  mol.GetCoords( target_atom , t_cds );
  mol.GetCoords( base_atom , b_cds );
  mol.GetCoords( x_atom , x_cds );

  float axis[3] , base_vec[3];
  DACLIB::join_vector( t_cds , b_cds , axis );
  DACLIB::join_vector( x_cds , b_cds , base_vec );
  float cdir[3];
  cdir[0] = dir[0]; cdir[1] = dir[1]; cdir[2] = dir[2];
  double a = -DACLIB::angle_about_axis( axis , base_vec , cdir );
  DACLIB::rotate_about_axis( axis , a , cdir );

  double ddir[3] = { cdir[0] , cdir[1] , cdir[2] };
  site.set_direction( ddir , GtplDefs::H_VECTOR , dir_num );

}

// ************************************************************************
void add_h_vector_to_o_site( OEMolBase &mol , OEAtomBase *target_atom ,
                             SinglePPhoreSite &site ) {

  if( count_exp_h_atoms_at_dist( mol , target_atom ) > 0 ) {

    // add the h atom using the explicit H - it seems churlish not to
    add_h_vector_using_exp_h( mol , target_atom , 1 , site );

  } else {

    OEAtomBase *base_atom;
    if( delocalised_o_or_n( target_atom , base_atom ) ) {
      add_planar_h_vector_to_o_site( mol , target_atom , base_atom , site );
    } else {
      add_tet_h_vector_to_o_site( mol , target_atom , site );
    }

  }

}

// ************************************************************************
void add_1_tet_h_vector_to_n_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                   SinglePPhoreSite &site ) {

  float att_vec[3];
  if( !get_attachment_vector( target_atom , mol , att_vec ) ) {
    return;
  }

  float td[15];
  make_oriented_tetrahedron( att_vec , td );

  // this one's a bit harder, as there's another attached atom that the first
  // direction needs to line up with, so that we get a proper tetrahedral
  // symmetry
  OEIter<OEAtomBase> atom = target_atom->GetAtoms( OEIsHeavy() );
  OEAtomBase *x = atom;
  if( !x ){
    return;
  }

  ++atom;
  OEAtomBase *xx = atom ? atom : (OEAtomBase *) 0;
  if( !xx ) {
    return;
  }

  // use a Diamond overlay to put the centre and first 2 pegs of the
  // tetrahedron on the target atom and 2 connected atoms, stick the
  // h atoms up the other leg.
  vector<float> at_cds( 9 , 0.0F );
  mol.GetCoords( target_atom , &at_cds[0] );
  mol.GetCoords( x , &at_cds[0] + 3 );
  mol.GetCoords( xx , &at_cds[0] + 6 );

  vector<float> td_cds( td , td + 9 );
  float l1 = DACLIB::distance( &at_cds[0] , &at_cds[3] );
  float l2 = DACLIB::distance( &at_cds[0] , &at_cds[6] );

  td_cds[3] *= l1; td_cds[4] *= l1; td_cds[5] *= l1;
  td_cds[6] *= l2; td_cds[7] *= l2; td_cds[8] *= l2;

  OverlayTrans ot( &td_cds[0] , &at_cds[0] , 3 );
  ot.overlay( 5 , td );

  // if the N atom is in a 5- or 6-membered ring, we want an equatorial N-H,
  // not an axial one. Do this by finding the in-ring atom attached to xx and
  // picking the vector from the two possible that has the lower dot-product.
  if( OEAtomIsInRingSize( target_atom , 5 ) ||
      OEAtomIsInRingSize( target_atom , 6 ) ||
      OEAtomIsInRingSize( target_atom , 7 ) ) {
    OEIter<OEAtomBase> xxx = xx->GetAtoms( OEAtomIsInRing() );
    if( xxx == target_atom ) {
      ++xxx;
    }
    float xx_cds[3] , xxx_cds[3] , xxx_vec[3] , v1[3] , v2[3];
    mol.GetCoords( xx , xx_cds );
    mol.GetCoords( xxx , xxx_cds );
    DACLIB::join_vector( xx_cds , xxx_cds , xxx_vec );
    DACLIB::join_vector( td , td + 9 , v1 );
    DACLIB::join_vector( td , td + 12 , v2 );
    DACLIB::normalise( xxx_vec );
    DACLIB::normalise( v1 );
    DACLIB::normalise( v2 );
    if( DACLIB::dot_product( xxx_vec , v1 ) <
        DACLIB::dot_product( xxx_vec , v2 ) ) {
      float d[3] = { td[9] - td[0] , td[10] - td[1] , td[11] - td[2] };
      site.set_direction( d , GtplDefs::H_VECTOR );
    } else {
      float d[3] = { td[12] - td[0] , td[13] - td[1] , td[14] - td[2] };
      site.set_direction( d , GtplDefs::H_VECTOR );
    }
  } else {
    // we have no reason to prefer one over the other
    float d[3] = { td[12] - td[0] , td[13] - td[1] , td[14] - td[2] };
    site.set_direction( d , GtplDefs::H_VECTOR );
  }

}

// ************************************************************************
void add_1_planar_h_vector_to_n_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                      SinglePPhoreSite &site ) {

  // get the two attached non-H atoms
  OEIter<OEAtomBase> atom = target_atom->GetAtoms( OEIsHeavy() );
  if( !atom ) return;
  OEAtomBase *x1 = 0 , *x2 = 0;
  x1 = atom;

  ++atom;
  if( !atom ) return;
  x2 = atom;

  float t_cds[3] , x1_cds[3] , x2_cds[3] , v1[3] , v2[3] , d[3];
  mol.GetCoords( target_atom , t_cds );
  mol.GetCoords( x1 , x1_cds );
  mol.GetCoords( x2 , x2_cds );

  DACLIB::join_vector( t_cds , x1_cds , v1 );
  DACLIB::join_vector( t_cds , x2_cds , v2 );
  d[0] = -v1[0] - v2[0];
  d[1] = -v1[1] - v2[1];
  d[2] = -v1[2] - v2[2];
  DACLIB::normalise( d );

  site.set_direction( d , GtplDefs::H_VECTOR );

}

// ************************************************************************
void add_1_h_vector_to_n_site( OEMolBase &mol , OEAtomBase *target_atom ,
                               SinglePPhoreSite &site ) {

  OEAtomBase *twiddle_base;

  if( delocalised_o_or_n( target_atom , twiddle_base ) ) {
    add_1_planar_h_vector_to_n_site( mol , target_atom , site );
  } else{
    add_1_tet_h_vector_to_n_site( mol , target_atom , site );
  }
}


// ************************************************************************
void add_2_tet_h_vectors_to_n_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                    SinglePPhoreSite &site ) {

  float att_vec[3];
  if( !get_attachment_vector( target_atom , mol , att_vec ) ) {
    return;
  }

  float td[15];
  make_oriented_tetrahedron( att_vec , td );

  float d[3];
  d[0] = -td[3];
  d[1] = -td[4];
  d[2] = -td[5];
  site.set_direction( d , GtplDefs::H_VECTOR );

  d[0] = -td[6];
  d[1] = -td[7];
  d[2] = -td[8];
  site.set_direction( d , GtplDefs::H_VECTOR );

}

// ************************************************************************
void add_1_tet_h_vector_to_nplus_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                       SinglePPhoreSite &site ) {

  // get the target coords and 3 coords of the 3 attached heavy atoms
  float tcds[3] , attcds[9];

  mol.GetCoords( target_atom , tcds );
  OEIter<OEAtomBase> atom = target_atom->GetAtoms( OEIsHeavy() );
  mol.GetCoords( atom , attcds );
  ++atom;
  mol.GetCoords( atom , attcds + 3 );
  ++atom;
  mol.GetCoords( atom , attcds + 6 );
  ++atom;

  // vectors from attached atoms to target atom
  float vecs[9];
  DACLIB::join_vector( attcds , tcds , vecs );
  DACLIB::join_vector( attcds + 3 , tcds , vecs + 3 );
  DACLIB::join_vector( attcds + 6 , tcds , vecs + 6 );

  // take sum of vecs and normalise, and thats the direction we want
  float d[3] = { vecs[0] , vecs[1] , vecs[2] };
  d[0] += vecs[3] + vecs[6];
  d[1] += vecs[4] + vecs[7];
  d[2] += vecs[5] + vecs[8];

  DACLIB::normalise( d );

  site.set_direction( d , GtplDefs::H_VECTOR );

}

// ************************************************************************
void add_2_tet_h_vectors_to_nplus_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                        SinglePPhoreSite &site ) {

  float att_vec[3];
  if( !get_attachment_vector( target_atom , mol , att_vec ) ) {
    return;
  }

  float td[15];
  make_oriented_tetrahedron( att_vec , td );

  OEIter<OEAtomBase> atom = target_atom->GetAtoms( OEIsHeavy() );
  OEAtomBase *x = atom;
  if( !x ) {
    return;
  }
  ++atom;
  OEAtomBase *xx = atom ? atom : (OEAtomBase *) 0;
  if( !xx ) {
    return;
  }

  // use a Diamond overlay to put the centre and first 2 pegs of the
  // tetrahedron on the target atom and 2 connected atoms, stick the
  // h atoms up the other 2 legs.
  vector<float> at_cds( 9 , 0.0F );
  mol.GetCoords( target_atom , &at_cds[0] );
  mol.GetCoords( x , &at_cds[0] + 3 );
  mol.GetCoords( xx , &at_cds[0] + 6 );

  vector<float> td_cds( td , td + 9 );
  float l1 = DACLIB::distance( &at_cds[0] , &at_cds[3] );
  float l2 = DACLIB::distance( &at_cds[0] , &at_cds[6] );

  td_cds[3] *= l1; td_cds[4] *= l1; td_cds[5] *= l1;
  td_cds[6] *= l2; td_cds[7] *= l2; td_cds[8] *= l2;

  OverlayTrans ot( &td_cds[0] , &at_cds[0] , 3 );
  ot.overlay( 5 , td );

  {
    float d[3] = { td[9] - td[0] , td[10] - td[1] , td[11] - td[2] };
    site.set_direction( d , GtplDefs::H_VECTOR );
  }
  {
    float d[3] = { td[12] - td[0] , td[13] - td[1] , td[14] - td[2] };
    site.set_direction( d , GtplDefs::H_VECTOR );
  }

}

// ************************************************************************
void add_3_tet_h_vectors_to_nplus_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                        SinglePPhoreSite &site ) {

  float att_vec[3];
  if( !get_attachment_vector( target_atom , mol , att_vec ) ) {
    return;
  }

  float td[15];
  make_oriented_tetrahedron( att_vec , td );

  float d[3];
  d[0] = -td[3];
  d[1] = -td[4];
  d[2] = -td[5];
  site.set_direction( d , GtplDefs::H_VECTOR );

  d[0] = -td[6];
  d[1] = -td[7];
  d[2] = -td[8];
  site.set_direction( d , GtplDefs::H_VECTOR );

  d[0] = -td[9];
  d[1] = -td[10];
  d[2] = -td[11];
  site.set_direction( d , GtplDefs::H_VECTOR );

}

// ************************************************************************
void add_2_planar_h_vectors_to_n_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                       OEAtomBase *base_atom ,
                                       SinglePPhoreSite &site ) {

  // get coords for target atom, base atom and an atom connected to base atom
  OEIter<OEAtomBase> atom = base_atom->GetAtoms( OEIsHeavy() );
  if( atom->GetIdx() == target_atom->GetIdx() ) {
    ++atom;
  }
  if( !atom ) {
    return; // wierd, though, what?
  }
  OEAtomBase *x_atom = atom;

  float b_cds[3] , t_cds[3] , x_cds[3];
  mol.GetCoords( target_atom , t_cds );
  mol.GetCoords( base_atom , b_cds );
  mol.GetCoords( x_atom , x_cds );

  // first H is parallel to x_atom, base_atom vector
  float x_b_vec[3] , att_vec[3];
  DACLIB::join_vector( x_cds , b_cds , x_b_vec );
  DACLIB::normalise( x_b_vec );
  site.set_direction( x_b_vec , GtplDefs::H_VECTOR );

  // 2nd H is 180 round from that
  DACLIB::join_vector( t_cds , b_cds , att_vec );
  DACLIB::rotate_about_axis( att_vec , M_PI , x_b_vec );
  site.set_direction( x_b_vec , GtplDefs::H_VECTOR );

}

// ************************************************************************
void add_2_h_vectors_to_n_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                SinglePPhoreSite &site ) {

  OEAtomBase *twiddle_base;
  if( delocalised_o_or_n( target_atom , twiddle_base ) ) {
    add_2_planar_h_vectors_to_n_site( mol , target_atom , twiddle_base ,
                                      site );
  } else {
    add_2_tet_h_vectors_to_n_site( mol , target_atom , site );
  }

}

// ************************************************************************
void add_tet_lp_vectors_to_o_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                   SinglePPhoreSite &site ) {

  float att_vec[3];
  if( !get_attachment_vector( target_atom , mol , att_vec ) ) {
    return;
  }

  float td[15];
  make_oriented_tetrahedron( att_vec , td );

  if( 1 == target_atom->GetHvyDegree() ) {

    float d[3];
    d[0] = -td[6];
    d[1] = -td[7];
    d[2] = -td[8];
    site.set_direction( d , GtplDefs::LP_VECTOR );

    d[0] = -td[9];
    d[1] = -td[10];
    d[2] = -td[11];
    site.set_direction( d , GtplDefs::LP_VECTOR );

  } else if( 2 == target_atom->GetHvyDegree() ) {

    OEIter<OEAtomBase> atom = target_atom->GetAtoms( OEIsHeavy() );
    OEAtomBase *x = atom;
    if( !x )
      return;
    ++atom;
    OEAtomBase *xx = atom ? atom : (OEAtomBase *) 0;
    if( !xx ) {
      return;
    }

    // use a Diamond overlay to put the centre and first 2 pegs of the
    // tetrahedron on the target atom and 2 connected atoms, stick the
    // lone pairs up the other 2 legs.
    vector<float> at_cds( 9 , 0.0F );
    mol.GetCoords( target_atom , &at_cds[0] );
    mol.GetCoords( x , &at_cds[0] + 3 );
    mol.GetCoords( xx , &at_cds[0] + 6 );

    vector<float> td_cds( td , td + 9 );
    float l1 = DACLIB::distance( &at_cds[0] , &at_cds[3] );
    float l2 = DACLIB::distance( &at_cds[0] , &at_cds[6] );

    td_cds[3] *= l1; td_cds[4] *= l1; td_cds[5] *= l1;
    td_cds[6] *= l2; td_cds[7] *= l2; td_cds[8] *= l2;

    OverlayTrans ot( &td_cds[0] , &at_cds[0] , 3 );
    ot.overlay( 5 , td );

    {
      float d[3] = { td[9] - td[0] , td[10] - td[1] , td[11] - td[2] };
      site.set_direction( d , GtplDefs::LP_VECTOR );
    }
    {
      float d[3] = { td[12] - td[0] , td[13] - td[1] , td[14] - td[2] };
      site.set_direction( d , GtplDefs::LP_VECTOR );
    }

  }

}

// ************************************************************************
void add_planar_lp_vectors_to_o_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                      SinglePPhoreSite &site ) {

  // two lone pairs, in the plane of the atoms attached to the atom
  // attached to target_atom;
  OEIter<OEAtomBase> atom = target_atom->GetAtoms( OEIsHeavy() );
  if( !atom ) {
    return; // bit of a small molecule, though
  }
  OEAtomBase *base_atom = atom;
  atom = base_atom->GetAtoms( OEIsHeavy() );
  if( atom->GetIdx() == target_atom->GetIdx() ) {
    ++atom;
  }
  if( !atom ) {
    return;
  }
  OEAtomBase *x_atom = atom;

  float b_cds[3] , t_cds[3] , x_cds[3];
  mol.GetCoords( target_atom , t_cds );
  mol.GetCoords( base_atom , b_cds );
  mol.GetCoords( x_atom , x_cds );

  float axis[3] , base_vec[3];
  DACLIB::join_vector( t_cds , b_cds , axis );
  DACLIB::join_vector( x_cds , b_cds , base_vec );
  DACLIB::normalise( base_vec );
  site.set_direction( base_vec , GtplDefs::LP_VECTOR );

  DACLIB::rotate_about_axis( axis , M_PI , base_vec );
  site.set_direction( base_vec , GtplDefs::LP_VECTOR );

}

// ************************************************************************
void add_lp_vectors_to_o_site( OEMolBase &mol , OEAtomBase *target_atom ,
                               SinglePPhoreSite &site ) {

  if( 2 == target_atom->GetTotalHCount() + target_atom->GetHvyDegree() ) {
    add_tet_lp_vectors_to_o_site( mol , target_atom , site );
  } else {
    add_planar_lp_vectors_to_o_site( mol , target_atom , site );
  }

}

// ************************************************************************
void add_lp_vector_to_aromatic_n_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                       SinglePPhoreSite &site ) {

  // straight out from ring. Also does substituted enamines.

  // get the two attached non-H atoms
  OEIter<OEAtomBase> atom = target_atom->GetAtoms( OEIsHeavy() );
  if( !atom ) return;
  OEAtomBase *x1 = 0 , *x2 = 0;
  x1 = atom;

  ++atom;
  if( !atom ) return;
  x2 = atom;

  float t_cds[3] , x1_cds[3] , x2_cds[3] , v1[3] , v2[3] , d[3];
  mol.GetCoords( target_atom , t_cds );
  mol.GetCoords( x1 , x1_cds );
  mol.GetCoords( x2 , x2_cds );

  DACLIB::join_vector( t_cds , x1_cds , v1 );
  DACLIB::join_vector( t_cds , x2_cds , v2 );
  d[0] = -v1[0] - v2[0];
  d[1] = -v1[1] - v2[1];
  d[2] = -v1[2] - v2[2];
  DACLIB::normalise( d );

  site.set_direction( d , GtplDefs::LP_VECTOR );

}

// ************************************************************************
void add_lp_vector_to_planar_n_site( OEMolBase &mol , OEAtomBase *target_atom ,
                                     SinglePPhoreSite &site ) {

  // want the base atom to be doubly bonded to target_atom - there might
  // be another heavy atom as well
  OEIter<OEBondBase> bond = target_atom->GetBonds( OEHasOrder( 2 ) );
  if( !bond ) {
    return;
  }

  OEAtomBase *base_atom = bond->GetEnd();
  if( base_atom->GetIdx() == target_atom->GetIdx() ) {
    base_atom = bond->GetBgn();
  }

  OEIter<OEAtomBase> xx_atom;
  OEAtomBase *x_atom = 0;
  for( xx_atom = base_atom->GetAtoms( OEIsHeavy() ) ; xx_atom ; ++xx_atom ) {
    if( xx_atom != target_atom ) {
      x_atom = xx_atom;
      break;
    }
  }
  if( !x_atom ) {
    return;
  }

  float b_cds[3] , t_cds[3] , x_cds[3];
  mol.GetCoords( target_atom , t_cds );
  mol.GetCoords( base_atom , b_cds );
  mol.GetCoords( x_atom , x_cds );

  float axis[3] , base_vec[3];
  DACLIB::join_vector( t_cds , b_cds , axis );
  DACLIB::join_vector( x_cds , b_cds , base_vec );
  DACLIB::normalise( base_vec );

  if( 1 == target_atom->GetTotalHCount() ) {
    DACLIB::rotate_about_axis( axis , M_PI , base_vec );
    site.set_direction( base_vec , GtplDefs::LP_VECTOR );
  }

}

// ************************************************************************
void add_lp_vectors_to_1_or_2_tet_n_site( OEMolBase &mol ,
                                          OEAtomBase *target_atom ,
                                          SinglePPhoreSite &site ) {

  // this is all a bit fiddly, as we need to be avoiding an H_VECTORS
  // that the atom will also have, from the corresponding donor site
  float att_vec[3];
  if( !get_attachment_vector( target_atom , mol , att_vec ) ) {
    return;
  }

  float td[15];
  make_oriented_tetrahedron( att_vec , td );

  if( 2 == target_atom->GetTotalHCount() ) {
    // just take the last peg in the tetrahedron, which won't be H_VECTORS
    float d[3];
    d[0] = -td[9];
    d[1] = -td[10];
    d[2] = -td[11];
    site.set_direction( d , GtplDefs::LP_VECTOR );
  } else if( 1 == target_atom->GetTotalHCount() ) {

    // this is all copied from add_1_tet_h_vector_to_n_site
    OEIter<OEAtomBase> atom = target_atom->GetAtoms( OEIsHeavy() );
    OEAtomBase *x = atom;
    if( !x ) {
      return;
    }

    ++atom;
    OEAtomBase *xx = atom ? atom : (OEAtomBase *) 0;
    if( !xx ) {
      return;
    }

    // use a Diamond overlay to put the centre and first 2 pegs of the
    // tetrahedron on the target atom and 2 connected atoms, stick the
    // lp up the other leg.
    vector<float> at_cds( 9 , 0.0F );
    mol.GetCoords( target_atom , &at_cds[0] );
    mol.GetCoords( x , &at_cds[0] + 3 );
    mol.GetCoords( xx , &at_cds[0] + 6 );

    vector<float> td_cds( td , td + 9 );
    float l1 = DACLIB::distance( &at_cds[0] , &at_cds[3] );
    float l2 = DACLIB::distance( &at_cds[0] , &at_cds[6] );

    td_cds[3] *= l1; td_cds[4] *= l1; td_cds[5] *= l1;
    td_cds[6] *= l2; td_cds[7] *= l2; td_cds[8] *= l2;

    OverlayTrans ot( &td_cds[0] , &at_cds[0] , 3 );
    ot.overlay( 5 , td );

    float d[3] = { td[12] - td[0] , td[13] - td[1] , td[14] - td[2] };
    site.set_direction( d , GtplDefs::LP_VECTOR );

  }

}

// ************************************************************************
void add_lp_vectors_to_3_tet_n_site( OEMolBase &mol ,
                                     OEAtomBase *target_atom ,
                                     SinglePPhoreSite &site ) {

  // we know there are 3 connected atoms
  OEIter<OEAtomBase> atom = target_atom->GetAtoms( OEIsHeavy() );
  OEAtomBase *x1 = atom;
  OEAtomBase *x2 = ++atom;
  OEAtomBase *x3 = ++atom;

  float t_cds[3] , x1_cds[3] , x2_cds[3] , x3_cds[3];
  float t_x1[3] , t_x2[3] , t_x3[3] , d[3];
  mol.GetCoords( target_atom , t_cds );
  mol.GetCoords( x1 , x1_cds );
  mol.GetCoords( x2 , x2_cds );
  mol.GetCoords( x3 , x3_cds );

  DACLIB::join_vector( t_cds , x1_cds , t_x1 );
  DACLIB::join_vector( t_cds , x2_cds , t_x2 );
  DACLIB::join_vector( t_cds , x3_cds , t_x3 );

  DACLIB::normalise( t_x1 );
  DACLIB::normalise( t_x2 );
  DACLIB::normalise( t_x3 );

  d[0] = - t_x1[0] - t_x2[0] - t_x3[0];
  d[1] = - t_x1[1] - t_x2[1] - t_x3[1];
  d[2] = - t_x1[2] - t_x2[2] - t_x3[2];
  DACLIB::normalise( d );

  site.set_direction( d , GtplDefs::LP_VECTOR );

}

// ************************************************************************
void add_lp_vectors_to_tet_n_site( OEMolBase &mol ,
                                   OEAtomBase *target_atom ,
                                   SinglePPhoreSite &site ) {

  if( 3 == target_atom->GetHvyDegree() )
    add_lp_vectors_to_3_tet_n_site( mol , target_atom , site );
  else
    add_lp_vectors_to_1_or_2_tet_n_site( mol , target_atom , site );

}

// ************************************************************************
void add_lp_vectors_to_nitrile_n_site( OEMolBase &mol ,
                                       OEAtomBase *target_atom ,
                                       SinglePPhoreSite &site ) {

  OEIter<OEAtomBase> base_atom = target_atom->GetAtoms( OEIsHeavy() );
  if( !base_atom )
    return;

  float t_cds[3] , b_cds[3] , d[3];
  mol.GetCoords( target_atom , t_cds );
  mol.GetCoords( base_atom , b_cds );
  DACLIB::join_vector( b_cds , t_cds , d );
  DACLIB::normalise( d );

  site.set_direction( d , GtplDefs::LP_VECTOR );

}

// ************************************************************************
void add_lp_vectors_to_n_site( OEMolBase &mol , OEAtomBase *target_atom ,
                               SinglePPhoreSite &site ) {

  if( 2 == target_atom->GetHvyDegree() &&
      0 == target_atom->GetTotalHCount() ) {
    add_lp_vector_to_aromatic_n_site( mol , target_atom , site );
  } else if( 1 == target_atom->GetHvyDegree() &&
             1 == target_atom->GetTotalHCount() ) {
    add_lp_vector_to_planar_n_site( mol , target_atom , site );
  } else if( 3 == target_atom->GetHvyDegree() + target_atom->GetTotalHCount() ) {
    add_lp_vectors_to_tet_n_site( mol , target_atom , site );
  } else if( 1 == target_atom->GetHvyDegree() &&
             0 == target_atom->GetTotalHCount() ) {
    add_lp_vectors_to_nitrile_n_site( mol , target_atom , site );
  }

}

// ************************************************************************
void make_12_virtual_points( SinglePPhoreSite &site ,
                             GtplDefs::DIRS_TYPE dir_type ) {

  if( !site.get_num_dirs() ) {
    return; // it's a bit odd, though
  }

  // get the 1st lone pair direction
  int dir_of_int = 0;
  for( int i = 0 ; i < site.get_num_dirs() ; ++i ) {
    if( dir_type == site.direction_type( i ) ) {
      dir_of_int = i;
      break;
    }
  }
  double lp_dir[3] = { site.direction( dir_of_int )[0] ,
                       site.direction( dir_of_int )[1] ,
                       site.direction( dir_of_int )[2] };
  DACLIB::normalise( lp_dir );
  double orig_lp_dir[3] = { lp_dir[0] , lp_dir[1] , lp_dir[2] };
  static const double rot_step = DACLIB::degrees_to_radians( 30.0 );
  double vs[36];
  for( int i = 0 ; i < 12 ; ++i ) {
    // virtual site is at distance 3Angstrom from atom
    vs[3*i] = 3.0 * lp_dir[0];
    vs[3*i+1] = 3.0 * lp_dir[1];
    vs[3*i+2] = 3.0 * lp_dir[2];

    lp_dir[0] = orig_lp_dir[0];
    lp_dir[1] = orig_lp_dir[1];
    lp_dir[2] = orig_lp_dir[2];
    DACLIB::rotate_about_axis( site.twiddle_axis() ,
                               double( i + 1 ) * rot_step , lp_dir );
  }

  site.set_virt_sites( vs , 12 );

}

//*****************************************************************************
void make_5_virtual_points( SinglePPhoreSite &site ,
                            GtplDefs::DIRS_TYPE dir_type ) {

  double outer_dirs[6];
  int dir_num = 0;
  // there will be 2 directions of dir_type, that's already confirmed
  // before this is called
  for( int i = 0 ; i < site.get_num_dirs() ; ++i ) {
    if( dir_type == site.direction_type( i ) ) {
      outer_dirs[3*dir_num] = site.direction( i )[0];
      outer_dirs[3*dir_num+1] = site.direction( i )[1];
      outer_dirs[3*dir_num+2] = site.direction( i )[2];
      ++dir_num;
    }
  }
  DACLIB::normalise( outer_dirs );
  DACLIB::normalise( outer_dirs + 3 );

  // the rotation axis is the cross product of the outer dirs.
  double rot_axis[3];
  DACLIB::norm_cross_product( outer_dirs , outer_dirs + 3 , rot_axis );

  double tot_angle = DACLIB::angle_about_axis( rot_axis , outer_dirs , outer_dirs + 3 );
  // virtual site is at distance 3 Angstrom from atom
  for( int i = 0 ; i < 6 ; ++i ) {
    outer_dirs[i] *= 3.0;
  }

  double vs[36];
  vs[0] = outer_dirs[0];
  vs[1] = outer_dirs[1];
  vs[2] = outer_dirs[2];

  double rot_step = tot_angle / 4.0;

  for( int i = 1 ; i < 4 ; ++i ) {
    double this_dir[3] = { outer_dirs[0] , outer_dirs[1] , outer_dirs[2] };
    DACLIB::rotate_about_axis( rot_axis , double( i ) * rot_step , this_dir );
    vs[3*i] = this_dir[0];
    vs[3*i+1] = this_dir[1];
    vs[3*i+2] = this_dir[2];
  }

  vs[12] = outer_dirs[3];
  vs[13] = outer_dirs[4];
  vs[14] = outer_dirs[5];

  site.set_virt_sites( vs , 5 );

}

//*****************************************************************************
// make a virtual point along each of the dirs of dir_type
void make_dir_virtual_points( SinglePPhoreSite &site ,
                              GtplDefs::DIRS_TYPE dir_type ) {

  double dir[3] , vs[36];

  int j = 0;
  for( int i = 0 ; i < site.get_num_dirs() ; ++i ) {
    if( dir_type == site.direction_type( i ) ) {
      dir[0] = site.direction( i )[0];
      dir[1] = site.direction( i )[1];
      dir[2] = site.direction( i )[2];
      DACLIB::normalise( dir );
      vs[3*j] = 3.0 * dir[0];
      vs[3*j+1] = 3.0 * dir[1];
      vs[3*j+2] = 3.0 * dir[2];
      ++j;
    }
  }

  if( site.get_flippable() && 1 == site.get_num_dirs() ) {
    // put another virtual site in flipped direction
    // flip axis is cross product of dir and twiddle_axis
    double flip_axis[3];
    DACLIB::cross_product( dir , site.twiddle_axis() , flip_axis );
    DACLIB::normalise( flip_axis );
    double ang1 = DACLIB::angle_about_axis( flip_axis , dir ,
                                            site.twiddle_axis() );
    double ang_to_rot = -2.0 * ( M_PI - ang1 );
    DACLIB::rotate_about_axis( flip_axis , ang_to_rot , dir );
    vs[3] = 3.0 * dir[0];
    vs[4] = 3.0 * dir[1];
    vs[5] = 3.0 * dir[2];
    j = 2;
  }

  site.set_virt_sites( vs , j );

}

// ************************************************************************
void add_virtual_acceptor_sites( OEMolBase &mol , OEAtomBase *target_atom ,
                                 SinglePPhoreSite &site ) {

  if( site.get_twiddlable() ) {
    // if the site is twiddlable, we need 12 virtual points in a cone
    // defined by the direction
    make_12_virtual_points( site , GtplDefs::LP_VECTOR );
  } else {
    if( 2 == site.get_num_dirs( GtplDefs::LP_VECTOR) ) {
      make_5_virtual_points( site , GtplDefs::LP_VECTOR );
    } else {
      make_dir_virtual_points( site , GtplDefs::LP_VECTOR );
    }
  }

}

// ************************************************************************
void add_virtual_donor_sites( OEMolBase &mol , OEAtomBase *target_atom ,
                              SinglePPhoreSite &site ) {

  if( site.get_twiddlable() ) {
    // if the site is twiddlable, we need 12 virtual points in a cone
    // defined by the direction
    make_12_virtual_points( site , GtplDefs::H_VECTOR );
  } else {
    make_dir_virtual_points( site , GtplDefs::H_VECTOR );
  }

}

// ************************************************************************
// create one or more SinglePPhoreSites with a direction along site_atom->H bond
// making as many different ones as there are H atoms on the site atom unless
// there's free rotation about the bond joining the site atom to the rest
// of the molecule. Assume there's only 1 site atom.
void create_pphore_sites_with_dir( GtplDefs::DIRS_TYPE dirs_type ,
                                   double *at_cds , const int *site_atoms ,
                                   int num_site_atoms , int type_code ,
                                   const string &type_string ,
                                   const string &label , OEMolBase &mol ,
                                   GtplDefs::DIR_MOVES dir_moves ,
                                   OEAtomBase *twiddle_base ,
                                   vector<SinglePPhoreSite *> &pharm_sites ) {

  // make the centroid, that's the easy bit.
  create_pphore_site( at_cds , site_atoms , num_site_atoms , type_code ,
                      type_string , label , root_mol_name( mol.GetTitle() ) ,
                      pharm_sites );

  SinglePPhoreSite *init_site = pharm_sites.back();
  init_site->set_dir_moves( dir_moves );

  // work out the twiddle axis if appropriate
  double dir[3] , dat_cds2[3];
  if( dir_moves != GtplDefs::NONE ) {
    mol.GetCoords( twiddle_base , dat_cds2 );
    DACLIB::join_vector( init_site->coords() , dat_cds2 , dir );
    float len = DACLIB::length( dir );
    dir[0] /= len; dir[1] /= len; dir[2] /= len;
    init_site->set_twiddle_axis( dir );
  }

  OEIter<OEAtomBase> atom;
  OEAtomBase *target_atom = 0;
  for( atom = mol.GetAtoms() ; atom ; ++atom ) {
    if( static_cast<unsigned int>( site_atoms[0] ) == atom->GetIdx() ) {
      target_atom = atom;
      break;
    }
  }

  target_atom = mol.GetAtom( OEHasAtomIdx( site_atoms[0] ) );

  if( GtplDefs::H_VECTOR == dirs_type ) {
    // add as many directions of type H_VECTOR as are appropriate
    if( 2 == target_atom->GetTotalHCount() && 1 == target_atom->GetHvyDegree() ) {
      add_2_h_vectors_to_n_site( mol , target_atom , *init_site );
    } else if( 1 == target_atom->GetTotalHCount() &&
               2 == target_atom->GetHvyDegree() ) {
      add_1_h_vector_to_n_site( mol , target_atom , *init_site );
    } else if( 1 == target_atom->GetTotalHCount() &&
               1 == target_atom->GetHvyDegree() ) {
      add_h_vector_to_o_site( mol , target_atom , *init_site );
    } else if( 1 == target_atom->GetTotalHCount() &&
               3 == target_atom->GetHvyDegree() ) {
      add_1_tet_h_vector_to_nplus_site( mol , target_atom , *init_site );
    } else if( 2 == target_atom->GetTotalHCount() &&
               2 == target_atom->GetHvyDegree() ) {
      add_2_tet_h_vectors_to_nplus_site( mol , target_atom , *init_site );
    } else if( 3 == target_atom->GetTotalHCount() &&
               1 == target_atom->GetHvyDegree() ) {
      add_3_tet_h_vectors_to_nplus_site( mol , target_atom , *init_site );
    }
    add_virtual_donor_sites( mol , target_atom , *init_site );
  } else if( GtplDefs::LP_VECTOR == dirs_type ) {
    if( OEElemNo::O == target_atom->GetAtomicNum() ||
        OEElemNo::S == target_atom->GetAtomicNum() ) {
      add_lp_vectors_to_o_site( mol , target_atom , *init_site );
    } else if( OEElemNo::N == target_atom->GetAtomicNum() ||
               OEElemNo::P == target_atom->GetAtomicNum() ) {
      add_lp_vectors_to_n_site( mol , target_atom , *init_site );
    }
    add_virtual_acceptor_sites( mol , target_atom , *init_site );
  }

}

// ****************************************************************************
void do_topol_category( OEMolBase &mol , OESubSearch &cat , float cat_factor ,
                        vector<float> &topol_factors ) {

  OEIter<OEMatchBase> match , match_start;

  match_start = cat.Match( mol );
  for( match = match_start ; match ; ++match ) {
    // the SMARTS are all set up only to return 1 atom, which makes life a bit
    // easier
    OEIter<OEMatchPair<OEAtomBase> > mp = match->GetAtoms();
    int this_at = mp->target->GetIdx();
    topol_factors[this_at] *= cat_factor;
  }

}

// ****************************************************************************
void do_topol_categories9_to_12( OEMolBase &mol , OESubSearch &cat9 ,
                                 OESubSearch &cat10 , OESubSearch &cat11 ,
                                 float cat_factor9 , float cat_factor10 ,
                                 float cat_factor11 , float cat_factor12 ,
                                 vector<float> &topol_factors ) {

  const int max_idx = mol.GetMaxAtomIdx();
  vector<char> atoms_9( max_idx , 0 );
  vector<char> atoms_10( max_idx , 0 );
  vector<char> atoms_11( max_idx , 0 );

  OEIter<OEMatchBase> match , match_start;

  match_start = cat9.Match( mol );
  for( match = match_start ; match ; ++match ) {
    OEIter<OEMatchPair<OEAtomBase> > mp = match->GetAtoms();
    atoms_9[mp->target->GetIdx()]++;
  }

  match_start = cat10.Match( mol );
  for( match = match_start ; match ; ++match ) {
    OEIter<OEMatchPair<OEAtomBase> > mp = match->GetAtoms();
    atoms_10[mp->target->GetIdx()]++;
  }

  match_start = cat11.Match( mol );
  for( match = match_start ; match ; ++match ) {
    OEIter<OEMatchPair<OEAtomBase> > mp = match->GetAtoms();
    atoms_11[mp->target->GetIdx()]++;
  }

  // any atoms in more than 1 of atoms_9 , atoms_10, atoms_11 are in cat 12,
  // otherwise they're in their own category
  for( int i = 0 ; i < max_idx ; ++i ) {
    int sum = atoms_9[i] + atoms_10[i] + atoms_11[i];
    if( sum >= 2 )
      topol_factors[i] *= cat_factor12;
    else {
      if( atoms_9[i] ) topol_factors[i] *= cat_factor9;
      if( atoms_10[i] ) topol_factors[i] *= cat_factor10;
      if( atoms_11[i] ) topol_factors[i] *= cat_factor11;
    }
  }

}

// ****************************************************************************
void do_topol_categories13_and_14( OEMolBase &mol , OESubSearch &cat13 ,
                                   float cat_factor13 , float cat_factor14 ,
                                   vector<float> &topol_factors ) {

  const int max_idx = mol.GetMaxAtomIdx();
  vector<char> atoms_13( max_idx , 0 );
  OEIter<OEMatchBase> match , match_start;

  match_start = cat13.Match( mol , true );
  for( match = match_start ; match ; ++match ) {
    OEIter<OEMatchPair<OEAtomBase> > mp = match->GetAtoms();
    atoms_13[mp->target->GetIdx()]++;
  }

  for( int i = 0 ; i < max_idx ; ++i ) {
    if( atoms_13[i] == 1 )
      topol_factors[i] *= cat_factor13;
    else if( atoms_13[i] > 1 )
      topol_factors[i] *= cat_factor14;
  }

}

// ****************************************************************************
void find_hphobe_itmoc_topol_factors( OEMolBase &mol ,
                                      vector<float> &topol_factors ) {

  // first one is a dummy, as the categories count from 1
  static const float cat_factors[15] = { -1.0f , 0.0f , 0.0f , 0.0f , 0.0f ,
                                         0.0f , 0.0f , 0.0f , 0.0f , 0.6f ,
                                         0.6f , 0.6f , 0.0f , 0.25f , 0.0f };

  // SMARTS targets for the topology-dependent hydrophobicity factors
  static OESubSearch cat1( "[$([#1]),$([#7]),$([#8])]" );
  static OESubSearch cat2( "[SH]" );
  static OESubSearch cat3( "[$(*~[+,++,-,--]),$(*~*~[+,++,-,--])]" );
  // v3.0 of loob put atoms connected to amide NH not in category 4, which I
  // think they should be.
  /*
    SMARTS for cat4 in SMARTS file format. They're then expanded for speed.
    nh [#7;H] 1 0
    oh [#8;H] 1 0
    nhoh [$nh,$oh] 1 0
    okamideat [*&D2,D3] 1 0
    nnotsing [$([$nhoh]!-*)] 1 0
    namide [$([$nhoh]-[$okamideat]!-*)] 1 0
    delocno [$nnotsing,$namide] 1 0
    notdelocno [$nhoh&!$delocno] 1 0
    cat4 [$([!$nhoh]-[$notdelocno]),$([!$nhoh]-*-[$notdelocno])] 1 1
  */
  //  static OESubSearch cat4( "[$(*~[$([$([NH;!$(*~*=,#,:*)])]),$([$([OH;!$(*~*=,#,:*)])])]),$(*~*~[$([$([NH;!$(*~*=,#,:*)])]),$([$([OH;!$(*~*=,#,:*)])])])]" );
  static OESubSearch cat4( "[$([!$([$([#7;H]),$([#8;H])])]-[$([$([$([#7;H]),$([#8;H])])&!$([$([$([$([$([#7;H]),$([#8;H])])]!-*)]),$([$([$([$([#7;H]),$([#8;H])])]-[$([*&D2,D3])]!-*)])])])]),$([!$([$([#7;H]),$([#8;H])])]-*-[$([$([$([#7;H]),$([#8;H])])&!$([$([$([$([$([#7;H]),$([#8;H])])]!-*)]),$([$([$([$([#7;H]),$([#8;H])])]-[$([*&D2,D3])]!-*)])])])])]" );
  static OESubSearch cat5( "[$(*~[SH;!$(*~*=,#,:*)])]" );
  static OESubSearch cat6( "[$(*=O),$(*~*=O)]" );
  static OESubSearch cat7( "[$([S&X3,S&X4,S&X5,S&X6]),$(*~$([S&X3,S&X4,S&X5,S&X6]))]" );
  static OESubSearch cat8( "[$(S=*)]" );
  static OESubSearch cat9( "*~*~*=O" );
  static OESubSearch cat10( "*~*~$([S&X3,S&X4,S&X5,S&X6])" );
  static OESubSearch cat11( "*~$(S=*)" );
  /*
    no [#7,#8] 1 0
    okamideat [!P&D2,D3] 1 0
    nnotsing [$([$no]!-*)] 1 0
    namide [$([$no]-[$okamideat]!-*)] 1 0
    delocno [$nnotsing,$namide] 1 0
    notdelocno [$no&!$delocno] 1 0
    cat13 [$([!$no]-[$notdelocno])&!$([!$no]-[$delocno])] 1 0
  */
  static OESubSearch cat13( "[$([!$([#7,#8])]-[$([$([#7,#8])&!$([$([$([$([#7,#8])]!-*)]),$([$([$([#7,#8])]-[$([!P&D2,D3])]!-*)])])])])&!$([!$([#7,#8])]-[$([$([$([$([#7,#8])]!-*)]),$([$([$([#7,#8])]-[$([!P&D2,D3])]!-*)])])])]" );

  std::fill( topol_factors.begin() , topol_factors.end() , 1.0F );

  // categories 1 to 8 are straightforward, topol factors being multiplied
  // together. 9 to 14 require a bit of counting.
  do_topol_category( mol , cat1 , cat_factors[1] , topol_factors );
  do_topol_category( mol , cat2 , cat_factors[2] , topol_factors );
  do_topol_category( mol , cat3 , cat_factors[3] , topol_factors );
  do_topol_category( mol , cat4 , cat_factors[4] , topol_factors );
  do_topol_category( mol , cat5 , cat_factors[5] , topol_factors );
  do_topol_category( mol , cat6 , cat_factors[6] , topol_factors );
  do_topol_category( mol , cat7 , cat_factors[7] , topol_factors );
  do_topol_category( mol , cat8 , cat_factors[8] , topol_factors );

  do_topol_categories9_to_12( mol , cat9 , cat10 , cat11 , cat_factors[9] ,
      cat_factors[10] , cat_factors[11] ,
      cat_factors[12] , topol_factors );
  do_topol_categories13_and_14( mol , cat13 , cat_factors[13] , cat_factors[14] ,
      topol_factors );

}

// ****************************************************************************
// special radii used by ITMOC calc in the past. Don't know where they came from,
// but they're different from those returned by DACLIB::atomic_num_to_rad.
// Put in for backward compatibility.
float atomic_num_to_loob_rad( int atomic_num ) {

  switch( atomic_num ) {
  case 1 : return 1.1;
  case 3 : return 1.1;
  case 6 : return 1.55;
  case 7 : return 1.45;
  case 8 : return 1.35;
  case 9 : return 1.3;
  case 11 : return 1.1;
  case 13 : return 2.0;
  case 14 : return 2.1;
  case 15 : return 1.75;
  case 16 : return 1.7;
  case 17 : return 1.65;
  case 19 : return 1.3;
  case 20 : return 1.5;
  case 35 : return 1.8;
  case 53 : return 2.05;
  }
  return 1.5F; // default value

}

// ****************************************************************************
// calculate the accesible surface area of those atoms that have a non-zero
// topological factor.  Since the final hydrophobic contribution is the product
// of the two, there's no point doing it if the topol_factor is zero.
void calc_h_factors( OEMolBase &mol , const vector<float> &topol_factors ,
                     bool include_aromatics , vector<char> &ats_to_consider ,
                     vector<float> &h_factors ) {

  const int max_idx = mol.GetMaxAtomIdx();
  vector<float> at_rads( max_idx , 0.0F );
  vector<float> at_rads_sq( max_idx , 0.0F );
  vector<float> at_cds( max_idx * 3 , 0.0F );
  vector<char> h_atom( max_idx , 0 );

  std::fill( ats_to_consider.begin() , ats_to_consider.end() , 1 );

  OEIter<OEAtomBase> atom;
  int at_num = 0;
  for( atom = mol.GetAtoms() ; atom ; ++atom , ++at_num ) {
    at_rads[at_num] = atomic_num_to_loob_rad( atom->GetAtomicNum() ) + 1.5F;
    at_rads_sq[at_num] = DACLIB::square( at_rads[at_num] );
    mol.GetCoords( atom , &at_cds[0] + 3 * at_num );
    if( 1 == atom->GetAtomicNum() ) {
      h_atom[at_num] = 1;
    }
  }
  int num_ats = at_num;

  at_num = 0;
  for( atom = mol.GetAtoms() ; atom ; ++atom , ++at_num ) {
    if( !include_aromatics && atom->IsAromatic() ) {
      ats_to_consider[at_num] = 0;
      continue;
    }
    if( fabs( topol_factors[at_num] ) < 1.0e-5 ) {
      h_factors[at_num] = 0.0F;
      continue;
    }

    float rad = at_rads[at_num];
    float this_at_cds[3] , pt_cds[3];
    this_at_cds[0] = at_cds[3 * at_num];
    this_at_cds[1] = at_cds[3 * at_num + 1];
    this_at_cds[2] = at_cds[3 * at_num + 2];

    int num_exp_pts = DACLIB::NUM_SPHERE_POINTS;
    for( int i = 0 ; i < DACLIB::NUM_SPHERE_POINTS ; ++i ) {
      pt_cds[0] = this_at_cds[0] + rad * DACLIB::SPHERE_POINTS[i][0];
      pt_cds[1] = this_at_cds[1] + rad * DACLIB::SPHERE_POINTS[i][1];
      pt_cds[2] = this_at_cds[2] + rad * DACLIB::SPHERE_POINTS[i][2];
      float *p_cds = &at_cds[0];
      for( int j = 0 ; j < num_ats ; ++j , p_cds += 3 ) {
        if( j == at_num || h_atom[j] )
          continue;
        float dist_sq = DACLIB::square( p_cds[0] - pt_cds[0] ) +
            DACLIB::square( p_cds[1] - pt_cds[1] ) +
            DACLIB::square( p_cds[2] - pt_cds[2] );
        if( dist_sq < at_rads_sq[j] ) {
          num_exp_pts--; // this point is occluded
          break;
        }
      }
    }
    h_factors[at_num] = float( num_exp_pts ) * DACLIB::square( rad - 1.5F ) *
        DACLIB::SPHERE_POINT_AREA * topol_factors[at_num];
  }

}

// ****************************************************************************
// test b1 is that all substituents to the ring are on one side of the plane
// of the ring. So can't apply to aromatics, where all substituents are in the
// plane.
bool ring_test_b1( OEMolBase &mol , const vector<int> &ring_ats ,
                   const vector<float> ring_at_cds ,
                   vector<OEAtomBase *> &matched_ats ) {

  vector<char> in_ring( mol.GetMaxAtomIdx() , 0 );
  for( unsigned int i = 0 , is = ring_ats.size() ; i < is ; ++i ) {
    in_ring[ring_ats[i]] = 1;
  }

  vector<float> subst_cds;
  OEIter<OEAtomBase> conn_at;
  float sat_cds[3];
  for( unsigned int i = 0 , is = matched_ats.size() ; i < is ; ++i ) {
    for( conn_at = matched_ats[i]->GetAtoms() ; conn_at ; ++conn_at ) {
      if( 1 != conn_at->GetAtomicNum() &&
          !in_ring[conn_at->GetIdx()] ) {
        subst_cds.push_back( ring_at_cds[3 * i] );
        subst_cds.push_back( ring_at_cds[3 * i + 1] );
        subst_cds.push_back( ring_at_cds[3 * i + 2] );
        mol.GetCoords( conn_at , sat_cds );
        subst_cds.push_back( sat_cds[0] );
        subst_cds.push_back( sat_cds[1] );
        subst_cds.push_back( sat_cds[2] );
      }
    }
  }

  // now have sets of coords for all substituents to the ring
  if( subst_cds.size() < 6 ) {
    return true; // 0 or 1 subst, so bound to be okay.
  }
  float vec1[3] , vec2[3];
  DACLIB::join_vector( &subst_cds[0] , &subst_cds[0] + 3 , vec1 );
  for( int i = 6 , is = subst_cds.size() ; i < is ; i += 6 ) {
    DACLIB::join_vector( &subst_cds[0] + i , &subst_cds[0] + i + 3 , vec2 );
    if( DACLIB::dot_product( vec1 , vec2 ) < 0.0F ) {
      return false; // in opposite dirs
    }
  }

  return true;

}

// ****************************************************************************
// ring test b2 is that at least 2 neighbouring ring atoms have h > 0 and no
// substituent of more than 2 atoms. ring_ats should come in with the first
// and last atom the same, to make my life easier. The definition is a verbatim
// quote from the paper, but it's ambiguous.  Do the two neighbouring ring atoms
// have to have no substituent of more than 2 atoms, or is that for the whole
// ring? I've gone with the former, and excluded H atoms from the defn of 2
// atoms.
bool ring_test_b2( OEMolBase &mol , const vector<int> &ring_ats ,
                   const vector<float> &h_factors ,
                   vector<OEAtomBase *> &matched_ats ) {

  vector<char> in_ring( mol.GetMaxAtomIdx() , 0 );
  for( unsigned int i = 0 , is = ring_ats.size() ; i < is ; ++i )
    in_ring[ring_ats[i]] = 1;

  for( unsigned int i = 0 , is = ring_ats.size() - 1 ; i < is ; ++i ) {
    if( h_factors[ring_ats[i]] > 0.0 && h_factors[ring_ats[i+1]] > 0.0 ) {
      // F atoms are added into the ring before ring_test_b2 is called, but
      // we don't want to deal with them here.
      if( OEElemNo::F == matched_ats[i]->GetAtomicNum() ||
          OEElemNo::F == matched_ats[i+1]->GetAtomicNum() ) {
        continue;
      }
      // find the first atoms in the molecule
      bool subst1 = false , subst2 = false;
      OEIter<OEAtomBase> conn_at , next_at;
      //      cout << "Substs on " << ring_ats[i] + 1 << endl;
      for( conn_at = matched_ats[i]->GetAtoms() ; conn_at ; ++conn_at ) {
        if( in_ring[conn_at->GetIdx()] )
          continue;
        if( conn_at->GetHvyDegree() > 2 ) {
          subst1 = true;
	}
        if( 2 == conn_at->GetHvyDegree() ) {
          // need to look at the next atom as well, which musn't be more than
          // GetHvyDegree 1.
          for( next_at = conn_at->GetAtoms() ; next_at ; ++next_at ) {
            if( next_at->GetIdx() != conn_at->GetIdx() &&
                next_at->GetHvyDegree() > 1 ) {
              subst1 = true;
	    }
          }
        }
      }
      // same again for 2nd atom. Probably should use a function
      //      cout << "Substs on " << ring_ats[i+1] + 1 << endl;
      for( conn_at = matched_ats[i+1]->GetAtoms() ; conn_at ; ++conn_at ) {
        if( in_ring[conn_at->GetIdx()] ) {
          continue;
	}
        if( conn_at->GetHvyDegree() > 2 ) {
          subst2 = true;
	}
        if( 2 == conn_at->GetHvyDegree() ) {
          // need to look at the next atom as well, which musn't be more than
          // GetHvyDegree 1.
          for( next_at = conn_at->GetAtoms() ; next_at ; ++next_at ) {
            if( next_at->GetIdx() != conn_at->GetIdx() &&
                next_at->GetHvyDegree() > 1 )
              subst2 = true;
          }
        }
      }
      //      cout << "Subst1 : " << subst1 << " and Subst2 : " << subst2 << endl;
      if( !subst1 && !subst2 ) {
        return true; // both are small enough so ring is ok.
      }
    }
  }

  //  cout << "Fails ring test b2" << endl;
  return false;

}

// ****************************************************************************
void make_itmoc_site( const vector<int> atoms , const vector<float> &cds ,
                      const vector<float> h_factors , const string &label ,
                      int type_num , float h_sum , const string &mol_name ,
                      vector<SinglePPhoreSite *> &pharm_sites ) {

  double pt_cds[3] = { 0.0 , 0.0 , 0.0 } , dir[3] = { 0.0 , 0.0 , 0.0 };
  for( unsigned int i = 0 ; i < atoms.size() ; ++i ) {
    float wt = h_factors[atoms[i]] / h_sum;
    pt_cds[0] += wt * cds[3 * i];
    pt_cds[1] += wt * cds[3 * i + 1];
    pt_cds[2] += wt * cds[3 * i + 2];
  }

  ostringstream oss;
  oss << "Site_" << pharm_sites.size() + 1 << "_" << label;
  try{
    SinglePPhoreSite *ps =
        new SinglePPhoreSite( pt_cds , dir , type_num , label , oss.str() ,
                              false , atoms.size() , &atoms[0] , mol_name );
    pharm_sites.push_back( ps );
  } catch( int max_atoms ) {
    // need to split it up, really, but just print an error for now and the point
    // won't be made.
    cerr << "site " << label << " in molecule " << mol_name
         << " has too many atoms : got " << atoms.size()
         << " but can't have more than " << MAX_SITE_ATOMS << endl;
  }

}

// ****************************************************************************
// do the ring itmoc sites - all rings up to 7 atoms
void make_ring_itmoc_hphobe_sites( const string &label , int type_num ,
                                   const string &parent_mol_name ,
                                   OEMolBase &mol ,
                                   const vector<float> &h_factors ,
                                   vector<char> &ats_to_consider ,
                                   vector<SinglePPhoreSite *> &pharm_sites ) {

  static OESubSearch rs[3] = { OESubSearch( "*1~*~*~*~*~1" ) ,
                               OESubSearch( "*1~*~*~*~*~*~1" ) ,
                               OESubSearch( "*1~*~*~*~*~*~*~1" ) };

  OEIter<OEMatchBase> match , match_start;
  OEIter<OEMatchPair<OEAtomBase> > mp;
  vector<int> ring_ats , all_ring_ats;
  vector<float> ring_at_cds;
  vector<OEAtomBase *> matched_ats;
  float at_cds[3];
  bool aromatic_ring;

  for( int i = 0 ; i < 3 ; ++i ) {
    match_start = rs[i].Match( mol , true );
    for( match = match_start ; match ; ++match ) {
      // if it's an aromatic ring, we may not want it
      matched_ats.clear();
      unsigned int aro_count = 0 , to_use = 0;
      for( mp = match->GetAtoms() ; mp ; ++mp ) {
        if( mp->target->IsAromatic() ) {
          ++aro_count;
	}
        if( ats_to_consider[mp->target->GetIdx()] ) {
          ++to_use;
	}
        matched_ats.push_back( mp->target );
      }
      if( to_use != matched_ats.size() ) {
        continue; // not using any of these atoms.
      }
      if( aro_count != matched_ats.size() ) {
        aromatic_ring = false;
      } else {
        aromatic_ring = true;
      }
      ring_ats.clear();
      ring_at_cds.clear();
      float h_sum = 0.0F;
      OEIter<OEAtomBase> nbour;
      for( unsigned int j = 0 , js = matched_ats.size() ; j < js ; ++j ) {
        int seq_num = matched_ats[j]->GetIdx();
        ring_ats.push_back( seq_num );
        mol.GetCoords( matched_ats[j] , at_cds );
        ring_at_cds.push_back( at_cds[0] );
        ring_at_cds.push_back( at_cds[1] );
        ring_at_cds.push_back( at_cds[2] );
        h_sum += h_factors[seq_num];
        // it seems that Catalyst includes F atoms connected to aromatic
        // rings in the ring
        if( aromatic_ring ) {
          for( nbour = matched_ats[j]->GetAtoms() ; nbour ; ++nbour ) {
            if( 9 == nbour->GetAtomicNum() ) {
              int seq_num = nbour->GetIdx();
              ring_ats.push_back( seq_num );
              matched_ats.push_back( nbour );
              mol.GetCoords( nbour , at_cds );
              ring_at_cds.push_back( at_cds[0] );
              ring_at_cds.push_back( at_cds[1] );
              ring_at_cds.push_back( at_cds[2] );
              h_sum += h_factors[seq_num];
            }
          }
        }
      }

      // mark these atoms to be taken out of consideration
      all_ring_ats.insert( all_ring_ats.end() , ring_ats.begin() , ring_ats.end() );

      if( h_sum < H_MIN ) {
        continue; // not greasy enough
      }
      // ring_test_b1 can't apply to aromatic rings, I don't think
      if( !aromatic_ring && ring_test_b1( mol , ring_ats , ring_at_cds ,
                                          matched_ats ) ) {
        // this is ok for a hphobe point
        make_itmoc_site( ring_ats , ring_at_cds , h_factors , label , type_num ,
                         h_sum , parent_mol_name , pharm_sites );
        continue;
      }
      ring_ats.push_back( ring_ats.front() ); // make it loop round
      matched_ats.push_back( matched_ats.front() );
      if( ring_test_b2( mol , ring_ats , h_factors , matched_ats ) ) {
        // this is ok for a hphobe point
        ring_ats.pop_back(); // to make it sensible again
        matched_ats.pop_back();
        make_itmoc_site( ring_ats , ring_at_cds , h_factors , label , type_num ,
                         h_sum , parent_mol_name , pharm_sites );
      }
    }
  }

  // take these atoms out of consideration
  for( unsigned int j = 0 , js = all_ring_ats.size() ; j < js ; ++j ) {
    ats_to_consider[all_ring_ats[j]] = 0;
  }

}

// ****************************************************************************
// do all groups for atoms with 3 or more bonds. A group is each atom with 3 or
// more bonds and those of its neighbours that are not bonded to any other atom
// provided that the sum of all h values is at least H_MIN.
void make_blob_itmoc_hphobe_sites( const string &label , int type_num ,
                                   const string &parent_mol_name ,
                                   OEMolBase &mol ,
                                   const vector<float> &h_factors ,
                                   vector<char> &ats_to_consider ,
                                   vector<SinglePPhoreSite *> &pharm_sites ) {

  OEIter<OEAtomBase> atom , nbour;
  int at_num = 0;
  vector<int> group_ats;
  vector<float> group_at_cds;
  float at_cds[3];

  for( atom = mol.GetAtoms() ; atom ; ++atom , ++at_num ) {
    if( !ats_to_consider[at_num] ) {
      continue;
    }
    if( atom->GetHvyDegree() > 2 ) {
      group_ats.clear();
      group_at_cds.clear();
      group_ats.push_back( atom->GetIdx() );
      mol.GetCoords( atom , at_cds );
      group_at_cds.push_back( at_cds[0] );
      group_at_cds.push_back( at_cds[1] );
      group_at_cds.push_back( at_cds[2] );
      float h_sum = h_factors[group_ats.back()];
      for( nbour = atom->GetAtoms() ; nbour ; ++nbour ) {
        if( 1 == nbour->GetHvyDegree() ) {
          group_ats.push_back( nbour->GetIdx() );
          mol.GetCoords( nbour , at_cds );
          group_at_cds.push_back( at_cds[0] );
          group_at_cds.push_back( at_cds[1] );
          group_at_cds.push_back( at_cds[2] );
          h_sum += h_factors[group_ats.back()];
        }
      }
      // take these atoms out of consideration
      for( unsigned int i = 0 , is = group_ats.size() ; i < is ; ++i ) {
        ats_to_consider[group_ats[i]] = 0;
      }

      if( h_sum > H_MIN ) {
        make_itmoc_site( group_ats , group_at_cds , h_factors , label ,
                         type_num , h_sum , parent_mol_name , pharm_sites );
      }
    }
  }

}

// ****************************************************************************
// build a contiguous chain of atoms from start_at that is pulled from atoms
// still available and such that the final chain has h factors between H_MIN
// and 2 * H_MIN.
float build_atom_chain( OEAtomBase *start_at , OEMolBase &mol ,
                        const vector<float> &h_factors ,
                        vector<char> &ats_to_consider , vector<int> &group_ats ,
                        vector<float> &group_at_cds ) {
  
  const static float two_h_min = 2 * H_MIN;

  group_ats.clear();
  group_at_cds.clear();
  float at_cds[3];
  group_ats.push_back( start_at->GetIdx() );
  ats_to_consider[group_ats.back()] = 0;
  mol.GetCoords( start_at , at_cds );
  group_at_cds.push_back( at_cds[0] );
  group_at_cds.push_back( at_cds[1] );
  group_at_cds.push_back( at_cds[2] );

  float h_sum = h_factors[group_ats.back()];
  // breadth first search for connected atoms that are still available.
  list<OEAtomBase *> still_to_do;
  still_to_do.push_back( start_at );
  int term_at = -1;
  if( 1 == start_at->GetHvyDegree() ) {
    term_at = 0;
  }

  OEIter<OEAtomBase> nbour;
  while( !still_to_do.empty() ) {
    OEAtomBase *next_at = still_to_do.front();
    still_to_do.pop_front();
    for( nbour = next_at->GetAtoms() ; nbour ; ++nbour ) {
      int seq_num = nbour->GetIdx();
      if( ats_to_consider[seq_num] ) {
        if( h_sum + h_factors[seq_num] < two_h_min ) {
          ats_to_consider[seq_num] = 0;
          still_to_do.push_back( nbour );
          group_ats.push_back( seq_num );
          mol.GetCoords( nbour , at_cds );
          group_at_cds.push_back( at_cds[0] );
          group_at_cds.push_back( at_cds[1] );
          group_at_cds.push_back( at_cds[2] );
          h_sum += h_factors[seq_num];
          if( 1 == nbour->GetHvyDegree() )
            term_at = group_ats.size() - 1;
        }
      }
    }
  }

  if( h_sum < H_MIN ) {
    group_ats.clear();
    group_at_cds.clear();
  }

  // if there's a terminal atom, set all atom coords to its, as we want the
  // hphobe point to be on that atom.
  if( -1 != term_at ) {
    term_at *= 3;
    for( unsigned int i = 0 , is = group_ats.size() ; i < is ; ) {
      group_at_cds[i++] = group_at_cds[term_at];
      group_at_cds[i++] = group_at_cds[term_at + 1];
      group_at_cds[i++] = group_at_cds[term_at + 2];
    }
  }

  return h_sum;

}

// ****************************************************************************
// all chains of enough atoms to make a bit greater than H_MIN and less than
// 2 * H_MIN.
void make_chain_itmoc_hphobe_sites( const string &label , int type_num ,
                                    const string &parent_mol_name , OEMolBase &mol ,
                                    const vector<float> &h_factors ,
                                    vector<char> &ats_to_consider ,
                                    vector<SinglePPhoreSite *> &pharm_sites ) {

  // first, cut out atoms with h_factor of 0.0F
  for( int i = 0 , is = h_factors.size() ; i < is ; ++i ) {
    if( fabs( h_factors[i] ) < 1.0e-10 )
      ats_to_consider[i] = 0;
  }

  while( 1 ) {

    // find the next starting atom in a chain, which will be the most connected,
    // highest sequence number atom that has yet to be done.
    int at_num = 0 , most_conns = -1;
    OEAtomBase *most_conn_at = 0;
    OEIter<OEAtomBase> atom;
    for( atom = mol.GetAtoms() ; atom ; ++atom , ++at_num ) {
      int deg = atom->GetHvyDegree();
      if( ats_to_consider[at_num] && deg >= most_conns ) {
        most_conns = deg;
        most_conn_at = atom;
      }
    }
    if( -1 == most_conns ) {
      break;
    }

    vector<int> group_ats;
    vector<float> group_at_cds;
    // build a chain out from this atom that is within h_sum limits.
    float h_sum = build_atom_chain( most_conn_at , mol , h_factors ,
                                    ats_to_consider , group_ats , group_at_cds );
    if( !group_ats.empty() && h_sum > H_MIN ) {
      make_itmoc_site( group_ats , group_at_cds , h_factors , label ,
                       type_num , h_sum , parent_mol_name , pharm_sites );
    }

  }

}

// ****************************************************************************
// Use the Catalyst method to make hydrophobe sites, as described in
// JCICS, 34, 1297-1308 (1994).
void make_itmoc_hphobe_pphore_sites( const string &label , int type_num ,
                                     OEMol &mol ,
                                     vector<vector<SinglePPhoreSite *> > &pharm_sites ,
                                     bool include_aromatics ) {

  vector<float> topol_factors( mol.GetMaxAtomIdx() , 1.0F );
  vector<float> h_factors( mol.GetMaxAtomIdx() , 0.0F );
  vector<char> ats_to_consider( mol.GetMaxAtomIdx() , 1 );
  OEIter<OEConfBase> conf;

  // the topology factors are geometry independent, so do once per molecule
  find_hphobe_itmoc_topol_factors( mol , topol_factors );

  int conf_num = 0;
  for( conf = mol.GetConfs() ; conf ; ++conf , ++conf_num ) {
    // cout << "conformation " << conf_num << " title : " << conf->GetTitle() << endl;
    // calculate the h factors for this conformation, re-setting ats_to_consider
    // at the same time, since the work for the previous conformation will
    // have changed it.
    calc_h_factors( conf , topol_factors , include_aromatics , ats_to_consider ,
                    h_factors );
    make_ring_itmoc_hphobe_sites( label , type_num , mol.GetTitle() , conf , h_factors ,
                                  ats_to_consider , pharm_sites[conf_num] );
    // do sites for atoms of 3 or more bonds
    make_blob_itmoc_hphobe_sites( label , type_num , mol.GetTitle() , conf , h_factors ,
                                  ats_to_consider , pharm_sites[conf_num] );
    // do chains
    make_chain_itmoc_hphobe_sites( label , type_num , mol.GetTitle() , conf , h_factors ,
                                   ats_to_consider , pharm_sites[conf_num] );
  }

}

// ****************************************************************************
GtplDefs::DIRS_TYPE assess_point_dir_type( PharmPoint &pharm_points ,
                                           const string &p_type ) {

  vector<string>::iterator s;
  s = find( pharm_points.ring_normal_points().begin() ,
            pharm_points.ring_normal_points().end() , p_type );
  if( s != pharm_points.ring_normal_points().end() ) {
    return GtplDefs::RING_NORMAL;
  }

  // or h direction?
  vector<string>::const_iterator q =
      find( pharm_points.h_vector_points().begin() ,
            pharm_points.h_vector_points().end() , p_type );
  if( q != pharm_points.h_vector_points().end() ) {
    return GtplDefs::H_VECTOR;
  }

  // or lone pair
  q = find( pharm_points.lp_vector_points().begin() ,
            pharm_points.lp_vector_points().end() , p_type );
  if( q != pharm_points.lp_vector_points().end() ) {
    return GtplDefs::LP_VECTOR;
  }

  return GtplDefs::UNKNOWN;

}

// ****************************************************************************
void make_pphore_sites( OEMol &mol , PharmPoint &pharm_points ,
                        const vector<pair<string,string> > &input_smarts ,
                        vector<pair<string,string> > &smarts_sub_defn ,
                        vector<vector<SinglePPhoreSite *> > &pharm_sites ) {

  // so we only build each SMARTS object once. There'll be a bit of a problem
  // if the SMARTS strings change from call to call, if there's no change
  // in the names.
  static map<string , OESubSearch *> smts_objs;

  // set up an empty set of SinglePPhoreSite vectors
  int num_confs = mol.NumConfs() , natoms = mol.NumAtoms();
  pharm_sites = vector<vector<SinglePPhoreSite *> >( num_confs ,
                                                     vector<SinglePPhoreSite *>() );

  double *at_cds = new double[3 * natoms];
  int    *at_nums = new int[natoms];

  map<string,vector<string> > &points_defs = pharm_points.points_defs();
  map<string,vector<string> >::iterator p;
  vector<string>::iterator q;
  map<string , OESubSearch *>::iterator r;
  int type_num = 0;
  string mol_name = mol.GetTitle();

  int itmoc_type_num = get_itmoc_type_num( pharm_points );

  // loop over each point type in pharm_points
  for( p = points_defs.begin() ; p != points_defs.end() ; ++p , ++type_num ) {
    if( p->second.empty() ) {
      continue; // point defined by key word (e.g. ITMOC, ITMOC_ALO) not SMARTS.
    }

    for( q = p->second.begin() ; q != p->second.end() ; ++q ) {
      r = smts_objs.find( *q );
      if( r == smts_objs.end() ) {
        // throws a DACLIB::SMARTSDefnError exception if there's a problem
        build_smarts_object( *q , input_smarts , smarts_sub_defn , smts_objs );
        r = smts_objs.find( *q );
      }
      // does it need a normal, h_vector or lone pair
      GtplDefs::DIRS_TYPE dirs_type = assess_point_dir_type( pharm_points , p->first );

      // match this SMARTS into the molecule, return only unique atom sets
      OEIter<OEMatchBase> match , match_start = r->second->Match( mol , true );
      OEIter<OEConfBase> conf;
      int conf_num = 0;
      
      // if this Point has an h direction, assess whether it's twiddlable - only
      // want to do this once per molecule, not once per conformation. We're
      // assuming that only 1 atom will be hit by the SMARTS
      vector<GtplDefs::DIR_MOVES> twiddlable;
      vector<OEAtomBase *> twiddle_base;
      if( dirs_type == GtplDefs::H_VECTOR ||
          dirs_type == GtplDefs::LP_VECTOR ) {
        for( match = match_start ; match ; ++match ) {
          OEIter<OEMatchPair<OEAtomBase> > mp = match->GetAtoms();
          OEAtomBase *twiddle_base_atom = 0;
          if( is_atom_twiddlable( mol , mp->target->GetIdx() ,
                                  twiddle_base_atom ) ) {
            twiddlable.push_back( GtplDefs::FREE );
          } else if( is_atom_flippable( mol , mp->target->GetIdx() ,
                                        twiddle_base_atom ) ) {
            twiddlable.push_back( GtplDefs::FLIP );
          } else {
            twiddlable.push_back( GtplDefs::NONE );
          }
          twiddle_base.push_back( twiddle_base_atom );
        }
      }
      
      for( conf = mol.GetConfs() ; conf ; ++conf , ++conf_num ) {
        int match_num = 0;
        for( match = match_start ; match ; ++match , ++match_num ) {
          string site_name = string( "Site_" ) + boost::lexical_cast<string>( pharm_sites[conf_num].size() + 1 ) +
              string( "_" ) + p->first;
          OEIter<OEMatchPair<OEAtomBase> > mp;
          double *pat_cds = at_cds;
          int    *pat_nums = at_nums;
          int num_hit_atoms = 0;
          for( mp = match->GetAtoms() ; mp ; ++mp ) {
            conf->GetCoords( mp->target , pat_cds );
            *pat_nums = mp->target->GetIdx();
            pat_cds += 3;
            ++pat_nums;
            ++num_hit_atoms;
          }
          switch( dirs_type ) {
          case GtplDefs::RING_NORMAL :
            create_pphore_site_with_normal( at_cds , at_nums , num_hit_atoms ,
                                            type_num , p->first , site_name ,
                                            mol_name , pharm_sites[conf_num] );
            break;
          case GtplDefs::H_VECTOR : case GtplDefs::LP_VECTOR :
            create_pphore_sites_with_dir( dirs_type , at_cds , at_nums ,
                                          num_hit_atoms , type_num ,
                                          p->first , site_name ,
                                          conf , twiddlable[match_num] ,
                                          twiddle_base[match_num] ,
                                          pharm_sites[conf_num] );
            break;
          default :
            create_pphore_site( at_cds , at_nums , num_hit_atoms , type_num ,
                                p->first , site_name , mol_name ,
                                pharm_sites[conf_num] );
            break;
          }
        }

      }
    }
  }

  delete [] at_nums;
  delete [] at_cds;

  // do hydrophobes the Catalyst way if required - hphobes_itmoc means aromatic
  // and aliphatic, hphobes_itmoc_alo means aliphatic only. They're mutually
  // exclusive.
  if( pharm_points.hphobes_itmoc() ) {
    make_itmoc_hphobe_pphore_sites( pharm_points.itmoc_label() , itmoc_type_num ,
                                    mol , pharm_sites , true );
  }
  if( pharm_points.hphobes_itmoc_alo() ) {
    make_itmoc_hphobe_pphore_sites( pharm_points.itmoc_alo_label() ,
                                    itmoc_type_num , mol , pharm_sites , false );
  }

}

// ****************************************************************************
void build_pphore_sites( vector<OEMol *> &oemols , PharmPoint &pharm_points ,
                         const vector<pair<string,string> > &input_smarts ,
                         vector<pair<string,string> > &smarts_sub_defn ,
                         vector<vector<vector<SinglePPhoreSite *> > > &pharm_sites ) {

  vector<OEMol *>::iterator p;
  vector<vector<SinglePPhoreSite *> > these_sites;
  for( p = oemols.begin() ; p != oemols.end() ; ++p ) {
    make_pphore_sites( *(*p) , pharm_points , input_smarts , smarts_sub_defn ,
                       these_sites );
    pharm_sites.push_back( these_sites );
    these_sites.clear();
  }

}

// ************************************************************************
// if two sites are on top of each other, jitter one so both can be seen
void jitter_pphore_sites( vector<BasePPhoreSite *> &pharm_sites ) {

  if( pharm_sites.size() <= 1 )
    return;

  vector<BasePPhoreSite *>::iterator p , q , pe , qe;
  pe = pharm_sites.end() - 1;
  qe = pharm_sites.end();

  double new_coords[3];
  for( p = pharm_sites.begin(); p != pe ; ++p ) {
    for( q = p + 1 ; q != qe ; ++q ) {
      const double *q_cds = (*q)->coords();
      if( (*p)->square_distance( q_cds ) < 2.5e-3 ) {
        new_coords[0] = q_cds[0] + drand48() * 0.01;
        new_coords[1] = q_cds[1] + drand48() * 0.01;
        new_coords[2] = q_cds[2] + drand48() * 0.01;
        (*q)->set_coords( new_coords );
      }
    }
  }

}

} // end of namespace DACLIB
