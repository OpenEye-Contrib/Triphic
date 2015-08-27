//
// file PluralityHit.cc
// David Cosgrove
// AstraZeneca
// 18th September 2007
//

#include "stddefs.H"
#include "OverlayTrans.H"
#include "PluralityHit.H"
#include "SinglePPhoreSite.H"

#include <limits>

#include <oechem.h>
#include <mpi.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;

// in eponymous file
OEMolBase *get_given_oeconf( OEMol &mol , int conf_num ,
                             bool add_conf_num_to_name );
void overlay_oemolbase( OEMolBase &mol , const OverlayTrans &ot );
// in eponymous file
BasePPhoreSite *read_site_from_file( istream &is );

namespace DACLIB {
// in grid_volumes.cc
VolumeGrid *prepare_mol_grid( OEMolBase *mol );
// in the eponymous fie
float atomic_num_to_rad( int atomic_num );
// in send_oemols_via_mpi.cc
void send_oemol_via_mpi( OEMolBase *mol , int dest_rank );
void rec_oemol_via_mpi( int source_rank , OEMolBase *mol );
void mpi_send_string( const string &str , int dest_rank );
void mpi_rec_string( int source_rank , string &str );
}

// *************************************************************************
PluralityHit::PluralityHit( OEMol &target_mol , int conf_num ,
                            vector<SinglePPhoreSite *> &target_sites ,
                            vector<int> &hit_sites ,
                            boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                            boost::shared_ptr<DACLIB::VolumeGrid> &soft_exc_vol_grid ,
                            vector<BTV> &hard_exc_vols ,
                            const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                            vector<BTC> &query_coords ,
                            vector<BTD> &query_dists ,
                            const vector<float> &hit_angles ,
                            const vector<float> &hit_torsions ) :
  ov_rms_( numeric_limits<float>::max() ) ,
  prot_overlap_( 0.0F ) , soft_exc_vol_( 0.0F ) {

  vector<float> t_cds;
  for( int i = 0 , is = hit_sites.size() ; i < is ; ++i ) {
    ov_sites_.push_back( boost::shared_ptr<SinglePPhoreSite>( new SinglePPhoreSite( *target_sites[hit_sites[i]] ) ) );
    t_cds.insert( t_cds.end() , target_sites[hit_sites[i]]->coords() ,
        target_sites[hit_sites[i]]->coords() + 3 );
  }

#ifdef NOTYET
  cout << endl << "Ov sites : " << endl;
  for( int i = 0 , is = ov_sites_.size() ; i < is ; ++i )
    cout << ov_sites_[i]->label() << " ";
  cout << endl;
#endif

  make_overlay_conf( target_mol , conf_num );
  if( !query_coords.empty() ) {
    overlay_hit( query_coords , t_cds );
  }

  if( !hard_exc_vols.empty() && !test_hard_exc_vols( hard_exc_vols ) ) {
#ifdef NOTYET
    cout << "ov_conf_ reset, not a hit" << endl;
#endif
    ov_conf_.reset();
    return; // it's not a hit after all
  }

  if( !query_dists.empty() ) {
    fill_distances( query_dists );
  }

  hit_angles_ = hit_angles;
  hit_torsions_ = hit_torsions_;

  boost::shared_ptr<DACLIB::VolumeGrid> mol_vol_grid;
  if( protein_grid ) {
    prot_overlap_ = calc_overlap_vol( *protein_grid , mol_vol_grid );
  }
  if( soft_exc_vol_grid ) {
    soft_exc_vol_ = calc_overlap_vol( *soft_exc_vol_grid , mol_vol_grid );
  }

  for( int i = 0 , is = score_vol_grids.size() ; i < is ; ++i ) {
    grid_vol_scores_.push_back( calc_overlap_vol( *(score_vol_grids[i].second) ,
                                                  mol_vol_grid ) );
  }

}

// *************************************************************************
PluralityHit::~PluralityHit() {

}

// *************************************************************************
void PluralityHit::write_to_stream( ostream &os , const string &sep ) {

  os << ov_conf_->GetTitle() << sep << ov_rms_ << sep << prot_overlap_
     << sep << soft_exc_vol_;
  for( int i = 0 , is = hit_dists_.size() ; i < is ; ++i )
    os << sep << hit_dists_[i];
  for( int i = 0 , is = hit_angles_.size() ; i < is ; ++i )
    os << sep << hit_angles_[i];
  for( int i = 0 , is = hit_torsions_.size() ; i < is ; ++i )
    os << sep << hit_torsions_[i];

  for( int i = 0 , is = grid_vol_scores_.size() ; i < is ; ++i ) {
    os << sep << grid_vol_scores_[i];
  }

#ifdef NOTYET
  for( int i = 0 , is = ov_sites_.size() ; i < is ; ++i ) {
    os << sep << ov_sites_[i]->label() << "(";
    for( int j = 0 , js = ov_sites_[i]->num_site_atoms() ; j < js ; ++j ) {
      os << ov_sites_[i]->site_atoms()[j];
      if( j != js - 1 )
        os << sep;
    }
    os << ")";
  }
#endif

  os << endl;

}

// *************************************************************************
bool PluralityHit::has_same_sites( PluralityHit &ph ) const {

  if( !ov_conf_ ) {
    return true;
  }

  vector<char> site_atoms( ov_conf_->GetMaxAtomIdx() , 0 );
  for( int i = 0 , is = ov_sites_.size() ; i < is ; ++i ) {
    fill( site_atoms.begin() , site_atoms.end() , 0 );
    int this_site_atoms = 0;
    for( int j = 0 , js = ov_sites_[i]->num_site_atoms() ; j < js ; ++j ) {
      site_atoms[ov_sites_[i]->site_atoms()[j]] = 1;
      ++this_site_atoms;
    }
    bool i_match_found = false;
    for( int j = 0 , js = ph.ov_sites_.size() ; j < js ; ++j ) {
      if( ov_sites_[i]->get_type_code() == ph.ov_sites_[j]->get_type_code() ) {
        int that_site_atoms = 0;
        for( int k = 0 , ks = ph.ov_sites_[j]->num_site_atoms() ; k < ks ; ++k ) {
          if( site_atoms[ph.ov_sites_[j]->site_atoms()[k]] )
            ++that_site_atoms;
        }
        if( this_site_atoms == that_site_atoms ) {
          i_match_found = true;
          break;
        }
      }
    }
    if( !i_match_found ) {
      return false;
    }
  }

  return true;

}

// *************************************************************************
void PluralityHit::add_scores_to_ov_conf( const vector<string> &grid_vol_files ) {

  if( !ov_conf_ ) {
    return;
  }

  OESetSDData( *ov_conf_ , "Grappel RMS" , lexical_cast<string>( ov_rms_ ) );
  OESetSDData( *ov_conf_ , "Grappel Excluded_Vol" ,
               lexical_cast<string>( prot_overlap_ ) );
  OESetSDData( *ov_conf_ , "Grappel Soft_Excluded_Vol" ,
               lexical_cast<string>( soft_exc_vol_ ) );

  for( int i = 0 , is = grid_vol_files.size() ; i < is ; ++i ) {
    boost::filesystem::path vp( grid_vol_files[i] );
    OESetSDData( *ov_conf_ , "Grappel " + vp.filename().string() ,
                 lexical_cast<string>( grid_vol_scores_[i] ) );
  }

}

// *************************************************************************
void PluralityHit::send_via_mpi( int dest_rank ) {

  unsigned int num_to_send = ov_conf_ && *ov_conf_ ? 1 : 0;
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , dest_rank , 0 , MPI_COMM_WORLD );
  if( num_to_send ) {
    DACLIB::send_oemol_via_mpi( ov_conf_.get() , dest_rank );
  }

  num_to_send = ov_sites_.size();
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , dest_rank , 0 , MPI_COMM_WORLD );
  ostringstream oss;
  for( unsigned int i = 0 ; i < num_to_send ; ++i ) {
    ov_sites_[i]->write_to_stream( oss );
  }
  DACLIB::mpi_send_string( oss.str() , dest_rank );

  MPI_Send( &ov_rms_ , 1 , MPI_FLOAT , dest_rank , 0 , MPI_COMM_WORLD );
  MPI_Send( &prot_overlap_ , 1 , MPI_FLOAT , dest_rank , 0 , MPI_COMM_WORLD );
  MPI_Send( &soft_exc_vol_ , 1 , MPI_FLOAT , dest_rank , 0 , MPI_COMM_WORLD );

  num_to_send = hit_dists_.size();
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , dest_rank , 0 , MPI_COMM_WORLD );
  if( num_to_send ) {
    MPI_Send( &hit_dists_[0] , num_to_send , MPI_FLOAT , dest_rank , 0 , MPI_COMM_WORLD );
  }

  num_to_send = hit_angles_.size();
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , dest_rank , 0 , MPI_COMM_WORLD );
  if( num_to_send ) {
    MPI_Send( &hit_angles_[0] , num_to_send , MPI_FLOAT , dest_rank , 0 , MPI_COMM_WORLD );
  }

  num_to_send = hit_torsions_.size();
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , dest_rank , 0 , MPI_COMM_WORLD );
  if( num_to_send ) {
    MPI_Send( &hit_torsions_[0] , num_to_send , MPI_FLOAT , dest_rank , 0 , MPI_COMM_WORLD );
  }

  num_to_send = grid_vol_scores_.size();
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , dest_rank , 0 , MPI_COMM_WORLD );
  if( num_to_send ) {
    MPI_Send( &grid_vol_scores_[0] , num_to_send , MPI_FLOAT , dest_rank , 0 , MPI_COMM_WORLD );
  }

}

// *************************************************************************
void PluralityHit::receive_via_mpi( int source_rank ) {

  unsigned int num_to_rec;
  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  if( num_to_rec ) {
    ov_conf_ = boost::shared_ptr<OEMolBase>( OENewMolBase( OEMolBaseType::OEDefault ) );
    DACLIB::rec_oemol_via_mpi( source_rank , ov_conf_.get() );
  }

  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  string str;
  DACLIB::mpi_rec_string( source_rank , str );
  istringstream iss( str );
  for( unsigned int i = 0 ; i < num_to_rec ; ++i ) {
    // read_site_from_file returns a BasePPhoreSite *, but we know that it's
    // a SinglePPhoreSite that has been packed in.
    ov_sites_.push_back( boost::shared_ptr<SinglePPhoreSite>( reinterpret_cast<SinglePPhoreSite *>( read_site_from_file( iss ) ) ) );
  }

  MPI_Recv( &ov_rms_ , 1 , MPI_FLOAT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &prot_overlap_ , 1 , MPI_FLOAT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &soft_exc_vol_ , 1 , MPI_FLOAT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );

  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  if( num_to_rec ) {
    hit_dists_ = vector<float>( num_to_rec , -1.0F );
    MPI_Recv( &hit_dists_[0] , num_to_rec , MPI_FLOAT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  }

  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  if( num_to_rec ) {
    hit_angles_ = vector<float>( num_to_rec , -1.0F );
    MPI_Recv( &hit_angles_[0] , num_to_rec , MPI_FLOAT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  }

  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  if( num_to_rec ) {
    hit_torsions_ = vector<float>( num_to_rec , -1.0F );
    MPI_Recv( &hit_torsions_[0] , num_to_rec , MPI_FLOAT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  }

  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  if( num_to_rec ) {
    grid_vol_scores_ = vector<float>( num_to_rec , -1.0F );
    MPI_Recv( &grid_vol_scores_[0] , num_to_rec , MPI_FLOAT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  }

}

// *************************************************************************
// make the overlay conformation, adding the hit site names as tagged info
void PluralityHit::make_overlay_conf( OEMol &target_mol , int conf_num ) {

  ov_conf_ = boost::shared_ptr<OEMolBase>( get_given_oeconf( target_mol ,
                                                             conf_num ,
                                                             true ) );

  string site_names;
  for( unsigned int i = 0 , is = ov_sites_.size() ; i < is ; ++i )
    site_names += ov_sites_[i]->get_full_name() + " ";

  // chop the last space off
  site_names = site_names.substr( 0 , site_names.length() - 1 );
  
  OESetSDData( *ov_conf_ , "Target_PPhore_Site_Names" , site_names );

}

// *************************************************************************
// overlay ov_conf_ and ov_sites_ using the transformation that moves
// targ_cds onto query_coords.
void PluralityHit::overlay_hit( vector<BTC> &query_coords ,
                                vector<float> &targ_cds ) {

  vector<float> q_cds;
  q_cds.reserve( query_coords.size() * 3 );
  for( int i = 0 , is = query_coords.size() ; i < is ; ++i ) {
    q_cds.push_back( query_coords[i].get<1>() );
    q_cds.push_back( query_coords[i].get<2>() );
    q_cds.push_back( query_coords[i].get<3>() );
  }

  OverlayTrans ot( &targ_cds[0] , &q_cds[0] , query_coords.size() );
  overlay_oemolbase( *ov_conf_ , ot );
  ov_rms_ = 0.0F;
  for( int i = 0 , is = ov_sites_.size() ; i < is ; ++i ) {
    overlay( *ov_sites_[i] , ot ); // in SinglePPhoreSite.cc
  }

}

// *************************************************************************
// make the distances for sites in the query and
// score the hit by the RMS difference between the actual pharm distance
// matrix for the given conformation and the mean of the required one
void PluralityHit::fill_distances( vector<BTD> &query_dists ) {

  ov_rms_ = 0.0F;
  for( int i = 0 , is = query_dists.size() ; i < is ; ++i ) {
    hit_dists_.push_back( ov_sites_[query_dists[i].get<0>()]->distance( ov_sites_[query_dists[i].get<1>()]->coords() ) );
    float mean_dist = 0.5 * ( query_dists[i].get<2>() + query_dists[i].get<3>() );
    ov_rms_ += DACLIB::square( hit_dists_.back() - mean_dist );
  }

  ov_rms_ = sqrt( ov_rms_ / float( query_dists.size() ) );

}

// *************************************************************************
// calculate the overlap volume between ov_conf_ and the given grid. Creates
// mol_grid if it doesn't already exist. Returns 0.0 if it can't
// do any better.
float PluralityHit::calc_overlap_vol( DACLIB::VolumeGrid &vol_grid ,
                                      boost::shared_ptr<DACLIB::VolumeGrid> &mol_grid ) {

  if( !ov_conf_ )
    return 0.0F;

  if( !mol_grid )
    mol_grid = boost::shared_ptr<DACLIB::VolumeGrid>( DACLIB::prepare_mol_grid( ov_conf_.get() ) );

  float sol_vol , surf_vol;
  mol_grid->common_volume( vol_grid , sol_vol , surf_vol );

  return sol_vol;

}

// *************************************************************************
// if any atom in ov_conf_ is within touching distance (distance apart less than
// sum of 2 radii) of an excluded vol sphere, return false as the hit fails.
bool PluralityHit::test_hard_exc_vols( vector<BTV> &hard_exc_vols ) {

  vector<float> sphere_cds , sphere_rads;
  for( int i = 0 , is = hard_exc_vols.size() ; i < is ; ++i ) {
    sphere_cds.push_back( hard_exc_vols[i].get<0>() );
    sphere_cds.push_back( hard_exc_vols[i].get<1>() );
    sphere_cds.push_back( hard_exc_vols[i].get<2>() );
    sphere_rads.push_back( hard_exc_vols[i].get<3>() );
  }

  OEIter<OEAtomBase> atom;
  float at_cds[3];
  for( atom = ov_conf_->GetAtoms() ; atom ; ++atom ) {
    float at_rad = DACLIB::atomic_num_to_rad( atom->GetAtomicNum() );
    ov_conf_->GetCoords( atom , at_cds );
    for( int i = 0 , is = sphere_rads.size() ; i < is ; ++i ) {
      if( DACLIB::sq_distance( &sphere_cds[3 * i] , at_cds ) <
          DACLIB::square( sphere_rads[i] + at_rad ) ) {
        return false;
      }
    }
  }

  return true;

}

