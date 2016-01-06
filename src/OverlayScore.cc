//
// file OverlayScore.cc
// David Cosgrove
// AstraZeneca
// 21st June 2006
//
// This files is the implementation file for the class OverlayScore, which holds the
// details of a clique based overlay between 2 conformations of 2 molecules
// and the scores for the overlay.

#include <limits>
#include <numeric>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>

#include "DiamondOverlay.H"
#include "OverlayScore.H"
#include "OverlayTrans.H"
#include "ShapeTanimoto.H"
#include "SinglePPhoreSite.H"
#include "VolumeGrid.H"

#include <oeszybki.h>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESz;

namespace DACLIB {
void grid_volumes( vector<OEMolBase *> &mols , float &include_vol ,
                   float &total_vol , float &shape_tanimoto );
void grid_surface_volumes( vector<OEMolBase *> &mols ,
                           float &include_surface_vol );
VolumeGrid *prepare_multi_mol_grid( vector<OEMolBase *> &mols );
VolumeGrid *prepare_mol_grid( OEMolBase *mol );
}

// these 2 in eponymous files
OEMolBase *get_given_oeconf( OEMol &mol , int conf_num ,
                             bool add_conf_num_to_name );
void overlay_oemolbase( OEMolBase &mol , const OverlayTrans &ot );

GtplDefs::SCORE_METHOD OverlayScore::score_method_ = GtplDefs::RMS_AND_SIZE;
boost::shared_ptr<OverlayScore> OverlayScore::ref_ov_;

// ******************************************************************************
OverlayScore::OverlayScore() : fixed_mol_name_( string( "Query" ) ) ,
  fixed_conf_( -1 ) , moving_conf_( -1 ) ,
  rms_( -1.0F ) , num_sites_( 0 ) ,
  hphobe_score_( -1.0F ) , hbond_score_( -1.0F ) , vol_score_( -1.0F ) ,
  included_vol_( -1.0F ) , total_vol_( -1.0F ) ,
  grid_shape_tanimoto_( -1.0F ) ,
  surface_vol_( -1.0F ) ,
  protein_clash_( -numeric_limits<float>::max() ) ,
  mmff_nrg_( numeric_limits<float>::max() ) ,
  num_sims_( 0 ) , clip_score_( -1.0F ) , clique_tanimoto_( -1.0F ) ,
  ov_trans_( 0 ) {

}

// ******************************************************************************
OverlayScore::OverlayScore( const string &fix_mol_name ,
                            const string &mov_mol_name ,
                            int fix_conf_num , int mov_conf_num ,
                            const vector<int> &clique ,
                            const vector<BasePPhoreSite *> &fixed_sites ,
                            const vector<SinglePPhoreSite *> &fixed_score_sites ,
                            const vector<SinglePPhoreSite *> &moving_sites ,
                            OEMolBase &fixed_conf , OEMol &moving_mol ,
                            boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid ,
                            boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                            const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                            boost::shared_ptr<OESzybki> &szybki,
                            bool use_ring_norms_in_overlay ,
                            bool use_h_vectors_in_overlay ,
                            bool use_lps_in_overlay , bool do_overlay ) :
  fixed_mol_name_( fix_mol_name ) , moving_mol_name_( mov_mol_name ) ,
  fixed_conf_( fix_conf_num ) , moving_conf_( mov_conf_num ) ,
  num_sites_( clique.size() / 2 ) ,
  hphobe_score_( -1.0F ) , hbond_score_( -1.0F ) , vol_score_( -1.0F ) ,
  included_vol_( 0.0F ) , total_vol_( 0.0F ) ,
  grid_shape_tanimoto_( -numeric_limits<float>::max() ) ,
  gauss_shape_tanimoto_( -numeric_limits<float>::max() ) ,
  surface_vol_( 0.0F ) ,
  protein_clash_( -numeric_limits<float>::max() ) ,
  mmff_nrg_( numeric_limits<float>::max() ) ,
  num_sims_( 0 ) , clip_score_( -1.0F ) , clique_tanimoto_( -1.0F ) ,
  ov_trans_( 0 ) {

  if( fixed_mol_name_.empty() ) {
    fixed_mol_name_ = string( "Query" );
  }

  if( num_sites_ > 50 ) {
    throw string( "Clique too large." );
  }

  copy( clique.begin() , clique.end() , site_nums_ );

  overlay_mol_and_sites( fixed_sites , moving_sites , moving_mol ,
                         szybki , use_ring_norms_in_overlay ,
                         use_h_vectors_in_overlay , use_lps_in_overlay ,
                         do_overlay );

  if( GtplDefs::GRID_SHAPE_TANI == get_score_method() ||
      GtplDefs::SURFACE_OVLP_VOLUME == get_score_method() ||
      GtplDefs::PROTEIN_CLASH == get_score_method() ) {
    calc_volume_scores( ov_conf_.get() , fixed_solid_grid , protein_grid ,
                        score_vol_grids );
  } else if( GtplDefs::GAUSS_SHAPE_TANI == get_score_method() ) {
    calc_gauss_shape_tanimoto( fixed_conf , *ov_conf_ );
  } else if( GtplDefs::ROBINS_SCORE_PARETO == get_score_method() ) {
    calc_robins_scores( fixed_score_sites , fixed_solid_grid );
  } else if( GtplDefs::OVERALL_SCORE_PARETO == get_score_method() ) {
    calc_volume_scores( ov_conf_.get() , fixed_solid_grid , protein_grid ,
                        score_vol_grids );
    calc_gauss_shape_tanimoto( fixed_conf , *ov_conf_ );
    calc_robins_scores( fixed_score_sites , fixed_solid_grid );
  }

}

// *************************************************************************
OverlayScore::OverlayScore( const std::string &fix_mol_name ,
                            const std::string &mov_mol_name ) :
  fixed_mol_name_( fix_mol_name ) , moving_mol_name_( mov_mol_name ) {

}

// *************************************************************************
// this one builds itself from a string, and is used, at least initially,
// in parallel triphic to create an object from data sent back from a slave.
// The string has most probably been written by write_contents_to_string().
OverlayScore::OverlayScore( const string &str ) :
  fixed_mol_name_( string( "Query" ) ) ,
  fixed_conf_( -1 ) , moving_conf_( -1 ) ,
  rms_( -1.0F ) , num_sites_( 0 ) ,
  hphobe_score_( -1.0F ) , hbond_score_( -1.0F ) , vol_score_( -1.0F ) ,
  included_vol_( -1.0F ) , total_vol_( -1.0F ) ,
  grid_shape_tanimoto_( -1.0F ) ,
  surface_vol_( -1.0F ) ,
  protein_clash_( -numeric_limits<float>::max() ) ,
  mmff_nrg_( numeric_limits<float>::max() ) ,
  num_sims_( 0 ) , clip_score_( -1.0F ) , clique_tanimoto_( -1.0F ) {

  istringstream iss( str );
  string tmp;
  iss >>  tmp;
  if( tmp != "<OverlayScore>" ) {
    throw std::runtime_error( "Error parsing OverlayScore string - unexpected " + tmp );
  }
  iss >> tmp;
  if( tmp != "<Details>" ) {
    throw std::runtime_error( "Error parsing OverlayScore string - unexpected " + tmp );
  }
  getline( iss , tmp ); // consume \n from previous read
  getline( iss , moving_mol_name_ ); // allow spaces in mol names
  iss >> moving_conf_;
  getline( iss , tmp ); // consume \n from previous read
  getline( iss , fixed_mol_name_ ); // allow spaces in mol names
  iss >> fixed_conf_;

  // there's a variable number of grid vols, from 0 upwards, at the end of this
  // line
  iss >> num_sites_ >> rms_
      >> hphobe_score_ >> hbond_score_ >> vol_score_ >> included_vol_
      >> total_vol_ >> grid_shape_tanimoto_ >> gauss_shape_tanimoto_
      >> surface_vol_ >> protein_clash_ >> clique_tanimoto_
      >> clip_score_ >> num_sims_ >> mmff_nrg_;

  int is = 0;
  iss >> is;
  for( int i = 0 ; i < is ; ++i ) {
    string grid_name;
    float grid_vol;
    iss >> grid_name >> grid_vol;
    grid_vols_.push_back( make_pair( grid_name , grid_vol ) );
  }

  for( int i = 0 ; i < 2 * num_sites_ ; ++i ) {
    iss >> site_nums_[i];
  }

  iss >> tmp;
  if( tmp != "</Details>" ) {
    throw std::runtime_error( "Error parsing OverlayScore string - unexpected " + tmp );
  }

  iss >> tmp;
  if( tmp == "<Molecule>" ) {

    oemolistream ims;
    size_t mol_start = str.find( "<Molecule>" );
    if( string::npos == mol_start )
      throw std::runtime_error( "Error parsing OverlayScore string - unexpected " + tmp );

    ims.SetFormat( OEFormat::SDF );
    ims.openstring( str.substr( mol_start + 11 ) );
    ov_conf_ =
      boost::shared_ptr<OEMolBase>( OENewMolBase( OEMolBaseType::OEDefault ) );
    ims >> *(ov_conf_.get());
    iss >> tmp; // should be the final </Molecule>

  }

  iss >> tmp; // should be the final </OverlayScore>

}

// ******************************************************************************
// this one builds itself from a string of scores from a scores output file,
// a conformation, assumed to be the overlaid one associated with the scores,
// a header line taken from a scores file and the separator string. It's
// for restarting a job from an intermediate output set.
OverlayScore::OverlayScore( boost::shared_ptr<OEMolBase> &hc , const string &scores_line ,
                            const string &headers_line , const string sep ,
                            bool no_hit_conf_number ) {

  ov_conf_ = hc; // the easy bit first

  vector<string> headers_splits , scores_splits;
  split_regex( headers_splits , headers_line , regex( sep ) );
  for( vector<string>::iterator p = headers_splits.begin() ; p != headers_splits.end() ; ++p ) {
    cout << *p << endl;
  }
  split_regex( scores_splits , scores_line , regex( sep ) );
  for( vector<string>::iterator p = scores_splits.begin() ; p != scores_splits.end() ; ++p ) {
    cout << *p << endl;
  }

  moving_mol_name_ = scores_splits[0];
  if( no_hit_conf_number ) {
    remove_conf_num_from_mol_name( moving_mol_name_ );
  }
  num_sites_ = lexical_cast<int>( scores_splits[1] );
  rms_ = lexical_cast<float>( scores_splits[2] );
  hphobe_score_ = lexical_cast<float>( scores_splits[3] );
  hbond_score_ = lexical_cast<float>( scores_splits[4] );
  vol_score_ = lexical_cast<float>( scores_splits[5] );
  included_vol_ = lexical_cast<float>( scores_splits[6] );
  total_vol_ = lexical_cast<float>( scores_splits[7] );
  grid_shape_tanimoto_ = lexical_cast<float>( scores_splits[8] );
  gauss_shape_tanimoto_ = lexical_cast<float>( scores_splits[9] );
  surface_vol_ = lexical_cast<float>( scores_splits[10] );
  protein_clash_ = lexical_cast<float>( scores_splits[11] );
  mmff_nrg_ = lexical_cast<float>( scores_splits[12] );
  clique_tanimoto_ = lexical_cast<float>( scores_splits[14] );
  clip_score_ = lexical_cast<float>( scores_splits[15] );
  num_sims_ = lexical_cast<int>( scores_splits[16] );
  int gvs = lexical_cast<int>( scores_splits[17] );
  grid_vols_ = vector<pair<string,float> >( gvs , make_pair( string( "" ) , 0.0F ) );
  for( int i = 0 ; i < gvs ; ++i ) {
    grid_vols_[i].first = headers_splits[18 + i];
    grid_vols_[i].second = lexical_cast<float>( scores_splits[18 + i] );
  }

}

// ******************************************************************************
OverlayScore::OverlayScore( const OverlayScore &ovs ) {

  copy_data( ovs );

}

// ******************************************************************************
OverlayScore::~OverlayScore() {

}

// ******************************************************************************
void OverlayScore::copy_data( const OverlayScore &ov ) {

  fixed_mol_name_ = ov.fixed_mol_name_;
  moving_mol_name_ = ov.moving_mol_name_;
  fixed_conf_ = ov.fixed_conf_;
  moving_conf_ = ov.moving_conf_;
  rms_ = ov.rms_;
  num_sites_ = ov.num_sites_;
  copy( ov.site_nums_ , ov.site_nums_ + 2 * num_sites_ , site_nums_ );
  hphobe_score_ = ov.hphobe_score_;
  hbond_score_ = ov.hbond_score_;
  vol_score_ = ov.vol_score_;
  included_vol_ = ov.included_vol_;
  total_vol_ = ov.total_vol_;
  grid_shape_tanimoto_ = ov.grid_shape_tanimoto_;
  gauss_shape_tanimoto_ = ov.gauss_shape_tanimoto_;
  surface_vol_ = ov.surface_vol_;
  num_sims_ = ov.num_sims_;
  protein_clash_ = ov.protein_clash_;
  mmff_nrg_ = ov.mmff_nrg_;
  clip_score_ = ov.clip_score_;
  clique_tanimoto_ = ov.clique_tanimoto_;
  grid_vols_ = ov.grid_vols_;

  if( ov.ov_trans_ ) {
    ov_trans_.reset( new OverlayTrans( *ov.ov_trans_ ) );
  } else {
    ov_trans_.reset();
  }

  ov_conf_ = ov.ov_conf_;
  ov_sites_ = ov.ov_sites_;

}

// ******************************************************************************
void OverlayScore::calc_volume_scores( boost::scoped_ptr<DACLIB::VolumeGrid> &mol_grid ,
                                       boost::shared_ptr<DACLIB::VolumeGrid> &target_solid_grid ,
                                       boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ) {

  if( protein_grid ) {
    calc_protein_clash( mol_grid.get() , protein_grid.get() );
  }

  if( !target_solid_grid ) {
    return; /* can't do anything else. This can happen if, for example, it's
         a Triphic sites-only query. */
  }

  mol_grid->common_volume( *target_solid_grid , included_vol_ , surface_vol_ );
  
  total_vol_ = mol_grid->solid_volume() + target_solid_grid->solid_volume() - included_vol_;
  grid_shape_tanimoto_ = included_vol_ / total_vol_;

}

// ****************************************************************************
// returns the final energy
void OverlayScore::optimise_overlay( boost::shared_ptr<OESzybki> &szybki ) {

  using namespace OESystem;
  OESzybkiResults sr;
  int natoms = ov_conf_->NumAtoms();
  OEIter<OEAtomBase> atom;

  bool imp_h = OEHasImplicitHydrogens( *ov_conf_ );
  vector<float> old_at_cds( 3 * natoms , 0.0F );
  int i = 0;
  for( atom = ov_conf_->GetAtoms( OEIsHeavy() ) ; atom ; ++atom , i +=3 ) {
    ov_conf_->GetCoords( atom , &old_at_cds[i] );
  }

  if( !(*szybki)( *ov_conf_ , sr ) ) {
    cerr << "Optimisation failed for " << ov_conf_->GetTitle() << endl;
    mmff_nrg_ = numeric_limits<float>::max();
  }

  // having moved the molecule, need to move the sites by the same
  // transformation
  vector<float> new_at_cds( 3 * natoms , 0.0F );
  i = 0;
  // Szybki adds explicit hydrogens, it seems, take them off if there weren't
  // there originally
  if( imp_h ) {
    OESuppressHydrogens( *ov_conf_ );
  }
  for( atom = ov_conf_->GetAtoms( OEIsHeavy() ) ; atom ; ++atom , i +=3 ) {
    ov_conf_->GetCoords( atom , &new_at_cds[i] );
  }
  OverlayTrans ov( &old_at_cds[0] , &new_at_cds[0] , natoms );

  overlay( ov_sites_ , ov );

  mmff_nrg_ = sr.GetTotalEnergy();

}

// ******************************************************************************
// make a copy of this OverlayScore, apart from the similars and the OverlayTrans
OverlayScore *OverlayScore::make_copy_no_sims() {

  OverlayScore *new_one = new OverlayScore( fixed_mol_name_ , moving_mol_name_ );
  new_one->copy_data( *this );
  new_one->num_sims_ = 0;

  return new_one;
  
}

// ******************************************************************************
// increment num_sims_.
void OverlayScore::add_similar() {

  ++num_sims_; // count this sim

}

// ******************************************************************************
// transfer the similars from sims to this.
void OverlayScore::add_similars( OverlayScore *sims ) {

  num_sims_ += sims->num_sims_;

}

// ******************************************************************************
void OverlayScore::clear_similars() {

  num_sims_ = 0;

}

// *************************************************************************
void OverlayScore::set_ov_conf( OEMolBase *new_ov_conf ) {

  ov_conf_ = boost::shared_ptr<OEMolBase>( new_ov_conf );

}

// *************************************************************************
OEMolBase *OverlayScore::get_ov_conf() {

  return ov_conf_.get();

}

// *************************************************************************
void OverlayScore::set_ov_sites( std::vector<SinglePPhoreSite *> &new_ov_sites ) {

  ov_sites_.clear();
  for( int i = 0 , is = new_ov_sites.size() ; i < is ; ++i ) {
    ov_sites_.push_back( boost::shared_ptr<SinglePPhoreSite>( new_ov_sites[i] ) );
  }

}

// ******************************************************************************
// calculate the volume based scores - using a grid-based method because it's
// easy. Gaussian volumes in the manner of ROCS etc. don't work well for multiple
// overlaid molecules.
void OverlayScore::calc_gauss_shape_tanimoto( OEMolBase &mol1 , OEMolBase &mol2 ) {

  DACLIB::ShapeTanimoto st( mol1 , mol2 );
  gauss_shape_tanimoto_ = st.shape_tanimoto();

}

// ******************************************************************************
void OverlayScore::calc_volume_scores( OEMolBase *mol ,
                                       boost::shared_ptr<DACLIB::VolumeGrid> &target_solid_grid ,
                                       boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ) {

  boost::scoped_ptr<DACLIB::VolumeGrid> g( DACLIB::prepare_mol_grid( mol ) );

  calc_volume_scores( g , target_solid_grid , protein_grid );

}

// ******************************************************************************
void OverlayScore::calc_volume_scores( OEMolBase *mol ,
                                       boost::shared_ptr<DACLIB::VolumeGrid> &target_solid_grid ,
                                       boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                                       const vector<pair<string,DACLIB::VolumeGrid *> > &vol_grids ) {

  boost::scoped_ptr<DACLIB::VolumeGrid> g( DACLIB::prepare_mol_grid( mol ) );

  calc_volume_scores( g , target_solid_grid , protein_grid );

  grid_vols_.clear();
  for( int i = 0 , is = vol_grids.size() ; i < is ; ++i ) {
    float sol_vol , surf_vol;
    vol_grids[i].second->common_volume( *g , sol_vol , surf_vol );
#ifdef NOTYET
    cout << i << " : " << vol_grids[i].first << " : " << sol_vol << " , " << surf_vol
         << " :: " << vol_grids[i].second->solid_volume() << " , "
         << vol_grids[i].second->surface_volume() << endl;
#endif
    grid_vols_.push_back( make_pair( vol_grids[i].first , sol_vol ) );
  }

}

// ******************************************************************************
// protein clash is 1.0 - the shape tanimoto between protein and mol, so giving
// higher values to lower overlaps. They'll be pretty close to 1.0, though, as
// the protein will be much bigger than the molecule.
void OverlayScore::calc_protein_clash( DACLIB::VolumeGrid *mol_grid ,
                                       DACLIB::VolumeGrid *prot_grid ) {

  float temp1 , temp2;
  if( prot_grid && mol_grid ) {
    mol_grid->common_volume( *prot_grid , temp1 , temp2 );
    protein_clash_ = temp1;
  }

}

// ******************************************************************************
void OverlayScore::overlay_mol_and_sites( const vector<BasePPhoreSite *> &fixed_sites ,
                                          const vector<SinglePPhoreSite *> &moving_sites ,
                                          OEMol &moving_mol ,
                                          boost::shared_ptr<OESzybki> &szybki,
                                          bool use_ring_norms_in_overlay ,
                                          bool use_h_vectors_in_overlay ,
                                          bool use_lps_in_overlay  ,
                                          bool do_overlay ) {

  vector<int> clique( site_nums_ , site_nums_ + 2 * num_sites_ );
  ov_trans_.reset( new OverlayTrans );
  try {
    rms_ = calc_overlay_trans( fixed_sites ,  moving_sites , clique , *ov_trans_ ,
                               false , false , false );
  } catch( DACLIB::DiamondOverlayError &e ) {
    throw OverlayScoreError( e.what() );
  }

  if( do_overlay ) {

    // the overlay process can move the directions on the fixed sites to work out
    // the best possible match of moving and fixed sites.  This is undesirable in a
    // parallel run as it can mean that different results are returned for the same
    // molecule depending on what other molecules that process has seen already. So
    // put them back at the end
    vector<GtplDefs::DIRS_TYPE> fixed_dir_types;
    vector<double> fixed_dirs( 9 * fixed_sites.size() , 0.0 );
    int dir_num = 0;
    for( int i = 0 , is = fixed_sites.size() ; i < is ; ++i ) {
      for( int j = 0 , js = fixed_sites[i]->get_num_dirs() ; j < js ; ++j , ++dir_num ) {
        fixed_dir_types.push_back( fixed_sites[i]->direction_type( j ) );
        fixed_dirs[3 * dir_num] = fixed_sites[i]->direction( j )[0];
        fixed_dirs[3 * dir_num + 1] = fixed_sites[i]->direction( j )[1];
        fixed_dirs[3 * dir_num + 2] = fixed_sites[i]->direction( j )[2];
      }
    }

    // overlay sites using ov_trans_ to make new set - moving_sites is not moved
    overlay_sites( fixed_sites , moving_sites , use_h_vectors_in_overlay ,
                   use_lps_in_overlay );

    // if needed, having moved the sites and lined up the h vectors and lone
    // pairs, calculate the transformation again to refine it using the moved
    // sites.
    if( use_h_vectors_in_overlay || use_lps_in_overlay ) {
      rms_ = calc_overlay_trans( fixed_sites ,  moving_sites , clique , *ov_trans_  ,
                                 use_ring_norms_in_overlay ,
                                 use_h_vectors_in_overlay , use_lps_in_overlay);
      overlay_sites( fixed_sites , moving_sites , use_h_vectors_in_overlay ,
                     use_lps_in_overlay );
    }
    // now overlay the relevant conformation of the moving_mol putting the answer
    // in ov_conf_, including info about the sites used
    overlay_moving_conf( moving_mol , fixed_sites , moving_sites );

    // optimise if that's required - don't do a single point calculation if
    // we're not doing a cutoff.
    if( szybki && OERunType::SinglePoint != szybki->GetRunType() ) {
      optimise_overlay( szybki );
      rms_ = calculate_site_rms( fixed_sites , ov_sites_ , clique );
    }

    // return directions to start values
    dir_num = 0;
    for( int i = 0 , is = fixed_sites.size() ; i < is ; ++i ) {
      for( int j = 0 , js = fixed_sites[i]->get_num_dirs() ; j < js ; ++j , ++dir_num ) {
        fixed_sites[i]->set_direction( &fixed_dirs[3 * dir_num] , fixed_dir_types[dir_num] , j );
      }
    }

  } else {
    ov_conf_ = boost::shared_ptr<OEMolBase>( get_given_oeconf( moving_mol ,
                                                        get_moving_conf() , false ) );
    for( int j = 0 , js = moving_sites.size() ; j < js ; ++j ) {
      ov_sites_.push_back( boost::shared_ptr<SinglePPhoreSite>( new SinglePPhoreSite( *moving_sites[j] ) ) );
    }
  }

}

// **************************************************************************
// free function that takes SinglePPhoreSites and overlays them using the
// OverlayTrans in the OverlayScore, returning copies of the sites, transformed
// to the overlay.
// Note that twiddle and flip can also change the direction vectors on the
// fixed sites, as they work out the best possible alignment of directions
// in the 2 overlays.
void OverlayScore::overlay_sites( const vector<BasePPhoreSite *> &fixed_sites ,
                                  const vector<SinglePPhoreSite *> &moving_sites ,
                                  bool use_h_vectors , bool use_lp_vectors ) {

  ov_sites_.clear();

  for( int j = 0 , js = moving_sites.size() ; j < js ; ++j ) {
    ov_sites_.push_back( boost::shared_ptr<SinglePPhoreSite>( new SinglePPhoreSite( *moving_sites[j] ) ) );
  }

  // in SinglePPhoreSite.cc
  overlay( ov_sites_ , *ov_trans() );

  const int *sits = sites();
  for( int i = 0 , is = 2 * num_sites() ; i < is ; i += 2 ) {
    if( use_h_vectors && ov_sites_[sits[i+1]]->get_twiddlable() ) {
      ov_sites_[sits[i+1]]->twiddle( *fixed_sites[sits[i]] ,
          GtplDefs::H_VECTOR );
    }
    if( use_lp_vectors && ov_sites_[sits[i+1]]->get_flippable() ) {
      ov_sites_[sits[i+1]]->flip( *fixed_sites[sits[i]] ,
          GtplDefs::LP_VECTOR );
    }
  }

}

// **************************************************************************
// free function that takes a multi-molecule query overlays it using the
// OverlayTrans in the OverlayScore, returning a new copy of target conf
// to the overlay.
void OverlayScore::overlay_moving_conf( OEMol &target_mol ,
                                        const vector<BasePPhoreSite *> &fixed_sites ,
                                        const vector<SinglePPhoreSite *> &moving_sites ) {

  ov_conf_ = boost::shared_ptr<OEMolBase>( get_given_oeconf( target_mol ,
                                                      get_moving_conf() , false ) );
  if( !ov_conf_ ) {
    return;
  }
  overlay_oemolbase( *ov_conf_ , *ov_trans() );

  add_clique_fixed_site_info( fixed_sites , *ov_conf_ );
  add_clique_moving_site_info( moving_sites , *ov_conf_ );

}

// ****************************************************************************
void OverlayScore::add_clique_fixed_site_info( const vector<BasePPhoreSite *> &fixed_sites ,
                                               OEMolBase &mol ) {

  string site_names;
  for( unsigned int i = 0 , is = num_sites() * 2 ; i < is ; i += 2 ) {
    site_names += fixed_sites[sites()[i]]->get_full_name() + " ";
  }

  // chop the last space off
  site_names = site_names.substr( 0 , site_names.length() - 1 );

  OESetSDData( mol , "Query_PPhore_Site_Names" , site_names );

}

// ****************************************************************************
void OverlayScore::add_clique_moving_site_info( const vector<SinglePPhoreSite *> &moving_sites ,
                                                OEMolBase &mol ) {

  string site_names;
  for( unsigned int i = 1 , is = num_sites() * 2 ; i < is ; i += 2 ) {
    site_names += moving_sites[sites()[i]]->get_full_name() + " ";
  }

  // chop the last space off
  site_names = site_names.substr( 0 , site_names.length() - 1 );

  OESetSDData( mol , "Target_PPhore_Site_Names" , site_names );

}

// *************************************************************************
// clique tanimoto is just the proportion of search sites that appear in the
// clique. Not sure anyone ever uses it, but triphic has always calculated it.
void OverlayScore::calc_clique_tanimoto( int num_query_sites ,
                                         int num_target_sites ) {

  float tot_sites( num_query_sites + num_target_sites );
  float f_clique_size( num_sites_ );

  clique_tanimoto_ = f_clique_size / ( tot_sites - f_clique_size );

}

// *************************************************************************
// calculate the clip score. Clip is a program similar to Triphic that Willet
// et al have published.  The clip score in this case is the clique tanimoto
// score times ( 1.0 - RMSD( matching distances ) / ( sum of max poss distance
// error )
void OverlayScore::calc_clip_score( vector<BasePPhoreSite *> &query_sites ,
                                    float dist_tol ) {

  if( clique_tanimoto_ < 0.0F ) {
    calc_clique_tanimoto( query_sites.size() , ov_sites_.size() );
  }

  float f_clique_size( num_sites_ );
  float a_am1 = f_clique_size * ( f_clique_size - 1.0F );
  float max_dist_err( a_am1 * dist_tol );

  float rmsd = 0.0F;
  for( int i = 0 ; i < num_sites_ - 1 ; ++i ) {
    for( int j = i + 1 ; j < num_sites_ ; ++j ) {
      float dist1 = DACLIB::distance( query_sites[site_nums_[2 * i]]->coords() ,
              query_sites[site_nums_[2 * j]]->coords() );
      float dist2 = DACLIB::distance( ov_sites_[site_nums_[2 * i + 1]]->coords() ,
              ov_sites_[site_nums_[2 * j + 1]]->coords() );
      rmsd += DACLIB::square( dist1 - dist2 );
    }
  }

  rmsd = sqrt( rmsd / a_am1 );

  clip_score_ = ( 1.0F - rmsd / max_dist_err ) * clique_tanimoto_;

}

// ******************************************************************************
// calculate the hphobe, donor and acceptor scores. The target volume grid is
// needed for calculating the occlusion of virtual donor and acceptor sites.
void OverlayScore::calc_robins_scores( const vector<SinglePPhoreSite *> &fixed_score_sites ,
                                       boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid ) {

  if( !fixed_solid_grid ) {
    return; // can't do Robin's scores as we need a volume for the occlusion part
  }

  vector<vector<SinglePPhoreSite *> > tmp_sites;
  tmp_sites.push_back( fixed_score_sites );
  tmp_sites.push_back( vector<SinglePPhoreSite *>() );
  transform( ov_sites_.begin() , ov_sites_.end() ,
             inserter( tmp_sites.back() , tmp_sites.back().begin() ) ,
             bind( &boost::shared_ptr<SinglePPhoreSite>::get , _1 ) );

  calc_directional_hphobe_score( tmp_sites );
  calc_hbond_score( tmp_sites , fixed_solid_grid );
  boost::shared_ptr<DACLIB::VolumeGrid> dummy_grid;
  calc_volume_scores( ov_conf_.get() , fixed_solid_grid , dummy_grid );
  vol_score_ = total_vol_;

#ifdef NOTYET
  cout << "Volume score = " << vol_score_ << endl;
#endif

}

// ******************************************************************************
// use the reference overlay to make a single number for Robin's Pareto ranking
float OverlayScore::calc_robins_pareto_score( const OverlayScore &os ) const {

  float ret_score = 0.0F;
  if( get_reference_overlay()->hphobe_score() != 0.0F ) {
#ifdef NOTYET
    cout << os.hphobe_score() << " / " << get_reference_overlay()->hphobe_score()
         << " = " << os.hphobe_score() / get_reference_overlay()->hphobe_score() << endl;
#endif
    ret_score += os.hphobe_score() / get_reference_overlay()->hphobe_score();
  }
  if( get_reference_overlay()->hbond_score() != 0.0F ) {
#ifdef NOTYET
    cout << os.hbond_score() << " / " << get_reference_overlay()->hbond_score()
         << " = " << os.hbond_score() / get_reference_overlay()->hbond_score() << endl;
#endif
    ret_score += os.hbond_score() / get_reference_overlay()->hbond_score();
  }
  // lower is better for volume score
  if( os.vol_score() != 0.0F ) {
#ifdef NOTYET
    cout << os.vol_score() << " / " << get_reference_overlay()->vol_score()
         << " = " <<  1.0F - os.vol_score() / get_reference_overlay()->vol_score() << endl;
    ret_score += 1.0F - os.vol_score() / get_reference_overlay()->vol_score();
    cout << "shape_tani = " << os.grid_shape_tanimoto() << endl;
#endif
    ret_score += os.grid_shape_tanimoto();
  }

  return ret_score;

}

// ******************************************************************************
// use the reference overlay to make a single number for Overall Pareto ranking
float OverlayScore::calc_overall_pareto_score( const OverlayScore &os ) const {

  float ret_score = 0.0F;

  // higher is better for num_sites()
  ret_score += float( os.num_sites() ) / float( get_reference_overlay()->num_sites() );
  // lower rms is better, the reference rms should be 0, so take 1.0 - rms().
  // It's possible that the rms() might be higher than 1.0, but it'll be in that
  // sort of ballpark, and a -ve overall score isn't a disaster.
  ret_score += 1.0 - os.rms();

  ret_score += calc_robins_pareto_score( os );

  // higher is better for included vol
  if( get_reference_overlay()->included_vol() != 0.0F && os.included_vol() != 0.0f ) {
    ret_score += os.included_vol() / get_reference_overlay()->included_vol();
  }

  // this is the same as vol_score so only count once. Lower is better.
  if( os.total_vol() != 0.0F && total_vol() != 0.0F  &&
      get_reference_overlay()->total_vol() && vol_score() < -0.5F ) {
    ret_score += 1.0F - os.total_vol() / get_reference_overlay()->total_vol();
  }

  // shape tanimotos are already normalised to 1.0, so just add
  if( os.grid_shape_tanimoto() > -numeric_limits<float>::max() ) {
    ret_score += os.grid_shape_tanimoto();
  }
  if( os.gauss_shape_tanimoto() > -numeric_limits<float>::max() ) {
    ret_score += os.gauss_shape_tanimoto();
  }

  // higher is better for surface volume
  if( get_reference_overlay()->surface_volume() != 0.0F && os.surface_volume() != 0.0F ) {
    ret_score += os.surface_volume() / get_reference_overlay()->surface_volume();
  }

  // the protein clash is in cubic angstrom which will swamp everything else.
  // Also, lower is better.  So doing 1.0 - fraction of molecule that clashes.
  if( get_reference_overlay()->protein_clash() > -numeric_limits<float>::max() ) {
    ret_score += 1.0 - ( os.protein_clash() / os.vol_score() );
  }

  // mmff_nrg is pretty much impossible to normalise, so ignore it.
  // clique tanimoto and clip score are also normalised to 1. They're set to
  // -1.0F if they haven't been calculated.
  if( os.clip_score() > -0.5F ) {
    ret_score += os.clip_score();
  }
  if( os.clique_tanimoto() > -0.5F ) {
    ret_score += os.clique_tanimoto();
  }

  return ret_score;

}

// ******************************************************************************
void OverlayScore::calc_directional_hphobe_score( const vector<vector<SinglePPhoreSite *> > &other_sites ) {

  vector<SiteCluster> clusters;
  // 1.5 is the distance threshold for the clustering.
  cluster_sites( "Hydrophobe" , 1.5 , other_sites , clusters );

  hphobe_score_ = score_hphobe_clusters( clusters );
#ifdef NOTYET
  cout << "Hphobe score : " << hphobe_score_ << endl;
#endif

}

// ******************************************************************************
void OverlayScore::calc_hbond_score( const vector<vector<SinglePPhoreSite *> > &other_sites ,
                                     boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid ) {

  vector<SiteCluster> donor_clusters;
  // 1.5 is the distance threshold for the clustering.
  cluster_sites( "Donor" , 1.5 , other_sites , donor_clusters );
  vector<float> donor_scores;
  score_hbond_clusters( donor_clusters , fixed_solid_grid , "Donor" , donor_scores );

  vector<SiteCluster> acc_clusters;
  // 1.5 is the distance threshold for the clustering.
  cluster_sites( "Acceptor" , 1.5 , other_sites , acc_clusters );
  vector<float> acc_scores;
  score_hbond_clusters( acc_clusters , fixed_solid_grid , "Acceptor" , acc_scores );

  // if a donor cluster has 75% of atoms in common with an acceptor cluster,
  // it's deemed a donor-acceptor cluster (probably hydroxyl groups) and only scored
  // once.
  remove_don_acc_clusters( donor_clusters , donor_scores ,
                           acc_clusters , acc_scores );

#ifdef NOTYET
  cout << "Number of clusters for scoring : " << donor_clusters.size() + acc_clusters.size()
       << endl;
#endif
  hbond_score_ = accumulate( donor_scores.begin() , donor_scores.end() , 0.0F ) +
      accumulate( acc_scores.begin() , acc_scores.end() , 0.0F );

#ifdef NOTYET
  cout << "HBond score = " << hbond_score_ << endl;
#endif

}

//*****************************************************************************
// cluster sites of given type not caring about whether their possession of a
// direction. Only 1 site in each cluster from each molecule
void OverlayScore::cluster_sites( const string &site_type , float thresh ,
                                  const vector<vector<SinglePPhoreSite *> > &other_sites ,
                                  vector<SiteCluster> &clusters ) {

  vector<SiteCluster> nnls;
  make_site_nnls( site_type , thresh , other_sites , nnls );

  cluster_nnls( nnls , clusters );

#ifdef NOTYET
  cout << "Clusters of sites of type " << site_type << endl;
  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    cout << "Cluster " << i << " : " << clusters[i].first << " :: ";
    for( int j = 0 , js = clusters[i].second.size() ; j < js ; ++j ) {
      cout << clusters[i].second[j].first << " , " << clusters[i].second[j].second->label() << " : ";
    }
    cout << endl;
  }
#endif

}

//*****************************************************************************
// make near-neighbour lists of sites with type string site_type
void OverlayScore::make_site_nnls( const string &site_type , float thresh ,
                                   const vector<vector<SinglePPhoreSite *> > &other_sites ,
                                   vector<SiteCluster> &nnls ) {

  // cout << "make_site_nnls : " << site_type << endl;

  float thresh_sq = thresh * thresh;

  // first generate all the near-neighbour lists (NNLs). For this, calculate the
  // distances from each site to all other sites, except those from the same
  // molecule, and keeping for each molecule the closest to the first site so
  // long as it's within the threshold
  for( int i = 0 , is = other_sites.size() ; i < is ; ++i ) {
    for( int j = 0 , js = other_sites[i].size() ; j < js ; ++j ) {
#ifdef NOTYET
      cout << other_sites[i][j]->get_type_string() << "XX" << site_type << "YY" << endl;
#endif
      if( !boost::iequals( other_sites[i][j]->get_type_string() , site_type ) ) {
       continue;
      }

      // cout << "nnl for " << other_sites[i][j]->label() << endl;
      vector<ClusterMember> next_nnl;
      next_nnl.push_back( make_pair( i , other_sites[i][j] ) );
      for( int k = 0 , ks = other_sites.size() ; k < ks ; ++k ) {
        if( k == i ) {
         continue;
        }
        float nearest_dist = numeric_limits<float>::max();
        int nearest_site = -1;
        for( int l = 0 , ls = other_sites[k].size() ; l < ls ; ++l ) {
          if( !boost::iequals( other_sites[k][l]->get_type_string() , site_type ) ) {
            continue;
          }
          float this_dist = other_sites[i][j]->square_distance( other_sites[k][l]->coords() );
          if( this_dist < nearest_dist ) {
            nearest_dist = this_dist;
            nearest_site = l;
          }
        }
        if( -1 != nearest_site && nearest_dist < thresh_sq ) {
          next_nnl.push_back( make_pair( k ,  other_sites[k][nearest_site] ) );
        }
      }
      nnls.push_back( make_pair( -1.0F , next_nnl ) );
    }
  }

#ifdef NOTYET
  cout << "Near neighbour lists for " << site_type << endl;
  for( int i = 0 , is = nnls.size() ; i < is ; ++i ) {
    for( int j = 0 , js = nnls[i].second.size() ; j < js ; ++j ) {
      cout << nnls[i].second[j].first << " , " << nnls[i].second[j].second->label() << " : ";
    }
    cout << endl;
  }
#endif

}

//*****************************************************************************
// calculate the sum of square distances of the cluster/nnl
void OverlayScore::calculate_nnl_dists( SiteCluster &nnl ) {

  float dist = 0.0F;
  for( int i = 0 , is = nnl.second.size() - 1 ; i < is ; ++i ) {
    for( int j = i + 1 , js = nnl.second.size() ; j < js ; ++j ) {
      dist += nnl.second[i].second->square_distance( nnl.second[j].second->coords() );
    }
  }

  nnl.first = dist;

}

//*****************************************************************************
// find the largest near-neighbour list, sum of square distances between members
// as a tie-breaker
int OverlayScore::find_largest_nnl( vector<SiteCluster> &nnls ) {

  int curr_largest = 0;
  for( int i = 1 , is = nnls.size() ; i < is ; ++i ) {
    if( nnls[i].second.size() > nnls[curr_largest].second.size() ) {
      curr_largest = i;
    } else if( nnls[i].second.size() == nnls[curr_largest].second.size() ) {
      if( nnls[i].first < 0.0F ) {
        calculate_nnl_dists( nnls[i] );
      }
      if( nnls[curr_largest].first < 0.0F ) {
        calculate_nnl_dists( nnls[curr_largest] );
      }
      if( nnls[i].first < nnls[curr_largest].first ) {
        curr_largest = i;
      }
    }
  }

  return curr_largest;

}
//*****************************************************************************
// put the near-neighbour lists into clusters. Destroys nnls in the process
void OverlayScore::cluster_nnls( vector<SiteCluster> &nnls ,
                                 vector<SiteCluster> &clusters )  {

  while( !nnls.empty() ) {
    //    cout << "nnls size : " << nnls.size() << endl;
    // next cluster is largest nnl, with distances as tie-breaker
    int next_nnl = find_largest_nnl( nnls );
    //    cout << "next_nnl : " << next_nnl << " : " << nnls[next_nnl].second.size() << endl;
    clusters.push_back( nnls[next_nnl] );
    // remove the members of nnls[next_nnl] from all other nnls
    for( int i = 0 , is = nnls.size() ; i < is ; ++i ) {
      if( i == next_nnl ) {
        continue;
      }
      for( int k = 0 , ks = nnls[next_nnl].second.size(); k < ks ; ++k ) {
        // if the seed of this nnl is in next_nnl, take out the whole nnl
        if( nnls[i].second[0] == nnls[next_nnl].second[k] ) {
          nnls[i].second.clear();
          continue;
        }
        for( int j = 1 , js = nnls[i].second.size() ; j < js ; ++j ) {
          if( nnls[i].second[j] == nnls[next_nnl].second[k] ) {
            nnls[i].second[j].second = static_cast<SinglePPhoreSite *>( 0 );
          }
        }
      }
    }
    // hose this nnl as it's captured
    nnls[next_nnl].second.clear();
    // close gaps
    for( int i = 0 , is = nnls.size() ; i < is ; ++i ) {
      int j = 0;
      for( int k = 0 , ks = nnls[i].second.size() ; k < ks ; ++k ) {
        if( nnls[i].second[k].second ) {
          nnls[i].second[j++] = nnls[i].second[k];
        }
      }
      nnls[i].second.erase( nnls[i].second.begin() + j , nnls[i].second.end() );
    }
    int j = 0;
    for( int i = 0 , is = nnls.size() ; i < is ; ++i ) {
      if( !nnls[i].second.empty() ) {
        nnls[j++] = nnls[i];
      }
    }
    nnls.erase( nnls.begin() + j , nnls.end() );
  }

}

//*****************************************************************************
// score the hpobe clusters, allowing for non-directional ones as well.
float OverlayScore::score_hphobe_clusters( const vector<SiteCluster> &clusters ) {

#ifdef NOTYET
  cout << "score_hphobe_clusters for " << clusters.size() << " clusters" << endl;
#endif

  float hy = 0.0F;
  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    // skip clusters of only 1 site, which the clustering algorithm can throw up
    if( 1 == clusters[i].second.size() ) {
#ifdef NOTYET
      cout << "skipping cluster " << i << " of size 1" << endl;
#endif
      continue;
    }
    // calculate the weighted mean square distance of each member from the centroid,
    // weight calculated according to Robin's formula (p454 of paper)
    float msd = mean_square_distance_to_centroid( clusters[i].second );
    float this_hy = 0.0F;
    if( msd <= 1.25 ) {
      this_hy += 1.0F - msd / 1.25F;
    }
    // and weighted mean of cosine between planes from same page. Only for those
    // sites that have directions
    float cp = mean_normal_cosine( clusters[i].second );
    float weighted_cp = cp >= 0.8 ? 10.0F * cp - 8.0F : 0.0;
    this_hy += weighted_cp;

#ifdef NOTYET
    cout << "cluster " << i << " size " << clusters[i].second.size()
         << " msd = " << msd << " weighted msd = " << 1.0F - msd / 1.25F << " cp = " << cp
         << " weighted cp = " << weighted_cp << " N2 = "
         << DACLIB::square( float( clusters[i].second.size() ) ) << " Score = "
         << DACLIB::square( float( clusters[i].second.size() ) ) * this_hy << endl;
#endif
    hy += DACLIB::square( float( clusters[i].second.size() ) ) * this_hy;
  }

  return hy;

}

//*****************************************************************************
void OverlayScore::score_hbond_clusters( vector<SiteCluster> &clusters ,
                                         boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid ,
                                         const string &hbond_type ,
                                         vector<float> &scores ) {

  scores.clear();

  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    // skip clusters of only 1 site, which the clustering algorithm can throw up
    if( 1 == clusters[i].second.size() ) {
      scores.push_back( -1.0F ); // to keep the clusters and score sizes the same
      continue;
    }
#ifdef NOTYET
    cout << "Scoring " << hbond_type << " cluster " << i << " of " << is << " size = " << clusters[i].second.size() << endl;
#endif
    vector<vector<SinglePPhoreSite *> > virt_don_sites;
    extract_virtual_sites( clusters[i].second , virt_don_sites );
    // cluster the virtual sites
    vector<SiteCluster> virt_clusters;
    cluster_sites( hbond_type , 1.5 , virt_don_sites , virt_clusters );
#ifdef NOTYET
    cout << "Number of virtual site clusters : " << virt_clusters.size() << endl;
#endif
    scores.push_back( calc_hbond_cluster_score( clusters[i] , virt_clusters , fixed_solid_grid ) );

    for( int j = 0 , js = virt_don_sites.size() ; j < js ; ++j ) {
      for( int k = 0 , ks = virt_don_sites[j].size() ; k < ks ; ++k ) {
        delete virt_don_sites[j][k];
      }
    }
  }

}

//*****************************************************************************
float OverlayScore::calc_hbond_cluster_score( SiteCluster &cluster ,
                                              vector<SiteCluster> &virt_clusters ,
                                              boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid ) {

  if( virt_clusters.empty() ) {
    // this might be the sign of something serious, but it might just be that
    // hbond acceptors and donors have been created without directions (the
    // default in older triphic versions) in which case there's nothing to be
    // done
    return 0.0F;
  }
  // a by-product of calc_occlusion is that the highest scoring (least occluded)
  // virtual sites cluster of those all of the same largest size, is put first
  // in virt_clusters
  float occ_score = calc_occlusion( cluster , virt_clusters , fixed_solid_grid );
  // The similarity factor sim is calculated from the strengths of the acceptor
  // groups in a cluster, which should be classified as strong, medium or weak.
  // This has yet to be implemented.
  float sim = 1.0F;
  if( cluster.first < 0.0F ) {
    calculate_nnl_dists( cluster );
  }
  float ap( cluster.second.size() );
  float fap = calc_fap( cluster.first );

  if( virt_clusters.front().first < 0.0F ) {
    calculate_nnl_dists( virt_clusters.front() );
  }
  float vp( virt_clusters.front().second.size() );
  float gvp = calc_gvp( virt_clusters.front().first );
  return( occ_score * sim * ( ap * ap * fap + vp * vp * gvp ) );

}

//*****************************************************************************
// calculate the occlusion score for the site cluster, using the clusters
// of virtual points.  If more than one cluster of virtual points has the same
// maximum size, do the score of them all and take the best. The best cluster
// is put at the top.
float OverlayScore::calc_occlusion( const SiteCluster &atom_sites ,
                                    vector<SiteCluster> &virt_sites ,
                                    boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid ) {

  // the occlusion score is 0.1 if all points below occluded, 1 if none are
  static const float occ_incr = 0.9 / 4;
  float occ = 0.0F;

  double atom_centroid[3];
  calc_cluster_centroid( atom_sites , atom_centroid );
#ifdef NOTYET
  cout << "Site centroid : " << atom_centroid[0] << " , " << atom_centroid[1]
       << " , " << atom_centroid[2] << endl;
#endif

  int best_clus = -1;
  for( int i = 0 , is = virt_sites.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "Virt sites cluster " << i << " of " << is << " size = " << virt_sites[i].second.size() << endl;
#endif
    if( virt_sites[i].second.size() < virt_sites.front().second.size() ) {
      break;
    }
    float this_occ = 1.0F;
    double virt_centroid[3];
    calc_cluster_centroid( virt_sites[i] , virt_centroid );
#ifdef NOTYET
    cout << "Virt sites centroid : " << virt_centroid[0] << " , "
         << virt_centroid[1] << " , " << virt_centroid[2] << endl;
#endif
    double vec[3];
    DACLIB::join_vector( atom_centroid , virt_centroid , vec );
    for( int j = 0 ; j < 4 ; ++j ) {
      float pos[3];
      pos[0] = atom_centroid[0] + vec[0] * ( 2.1 + float( j ) * 0.3 );
      pos[1] = atom_centroid[1] + vec[1] * ( 2.1 + float( j ) * 0.3 );
      pos[2] = atom_centroid[2] + vec[2] * ( 2.1 + float( j ) * 0.3 );
      if( fixed_solid_grid->point_in_volume( pos ) ) {
#ifdef NOTYET
        cout << "pos : " << pos[0] << " , " << pos[1] << " , " << pos[2] << " inside" << endl;
#endif
        this_occ -= occ_incr;
      }
    }
    if( this_occ > occ ) {
      occ = this_occ;
      best_clus = i;
    }
  }

#ifdef NOTYET
  cout << "num virt sites : " << virt_sites.front().second.size()
       << " occ score = " << occ << endl;
#endif
  if( best_clus > 0 ) {
    swap( virt_sites[0] , virt_sites[best_clus] );
  }
  return occ;

}

//*****************************************************************************
// calculate f(ap) as:
// fap falls linearly from 1.0 to 0.3 as ap increases from 0.15 to 0.75, OR
// 1.0 and 0.3 below and above these values
float OverlayScore::calc_fap( float ap ) {

  static const float m = -0.7F / 0.6F;
  if( ap < 0.15 ) {
    return 1.0;
  } else if( ap > 0.75 ) {
    return 0.3;
  } else {
    return ( m * ap ) + 1.175;
  }

}

//*****************************************************************************
// calculate g(vp) as:
// gvp falls linearly from 1.0 to 0.3 as vp increases from 0.5 to 1.5, OR
// 1.0 and 0.3 below and above these values
float OverlayScore::calc_gvp( float vp ) {

  static const float m = -0.7F;
  if( vp < 0.5 ) {
    return 1.0;
  } else if( vp > 1.5 ) {
    return 0.3;
  } else {
    return ( m * vp ) + 1.35;
  }

}

//*****************************************************************************
// if a donor cluster has 75% of atoms in common with an acceptor cluster,
// it's deemed a donor-acceptor cluster (probably hydroxyl groups) and only scored
// once so remove lower-scoring one.
void OverlayScore::remove_don_acc_clusters( std::vector<SiteCluster> &donor_clusters ,
                                            std::vector<float> &donor_scores ,
                                            std::vector<SiteCluster> &acc_clusters ,
                                            std::vector<float> &acc_scores ) const {

#ifdef NOTYET
  cout << "remove_don_acc_clusters : " << donor_clusters.size() << " , "
       << donor_scores.size() << " and " << acc_clusters.size()
       << " , " << acc_scores.size() << endl;
#endif

  for( unsigned int i = 0 ; i < donor_clusters.size() ; ++i ) {
    vector<ClusterMember> &dc = donor_clusters[i].second;
    for( unsigned int j = 0 ; j < acc_clusters.size() ; ++j ) {
      vector<ClusterMember> &ac = acc_clusters[j].second;
      int num_in_common = 0;
      for( int ii = 0 , iis = dc.size() ; ii < iis ; ++ii ) {
        for( int jj = 0 , jjs = ac.size() ; jj < jjs ; ++jj ) {
          // first is the molecule number, second is the site. See if the two sites
          // are in the same molecule and have the same site atom. They're donors
          // or acceptors, we know that, so there should only be 1 atom in each
          // site
          if( dc[ii].first == ac[jj].first &&
              dc[ii].second->site_atoms()[0] == ac[jj].second->site_atoms()[0] ) {
            ++num_in_common;
          }
        }
      }
      int num_at_75 = int( 0.75F * float( dc.size() + ac.size() ) );
#ifdef NOTYET
      cout << "remove_don_acc_clusters, num_in_common = " << num_in_common
           << " for donor cluster " << i << " acc cluster " << j << endl;
      cout << "75% of num atoms = " << num_at_75 << endl;
#endif
      if( 2 * num_in_common >= num_at_75 ) {
#ifdef NOTYET
        cout << "it's a common cluster : donor_score = " << donor_scores[i]
                << " and acc_score = " << acc_scores[j] << endl;
#endif
        if( donor_scores[i] < acc_scores[j] ) {
          dc.clear();
          donor_clusters[i].first = -1.0F;
          donor_scores[i] = -1.0F;
        } else {
          ac.clear();
          acc_clusters[j].first = -1.0F;
          acc_scores[j] = -1.0F;
        }
      }
    }
  }

  donor_clusters.erase( remove_if( donor_clusters.begin() , donor_clusters.end() ,
                                   bind( equal_to<float>() , -1.0F ,
                                         bind( &SiteCluster::first , _1 ) ) ) ,
                        donor_clusters.end() );
  donor_scores.erase( remove( donor_scores.begin() , donor_scores.end() , -1.0F ) ,
                      donor_scores.end() );
  acc_clusters.erase( remove_if( acc_clusters.begin() , acc_clusters.end() ,
                                 bind( equal_to<float>() , -1.0F ,
                                       bind( &SiteCluster::first , _1 ) ) ) ,
                      acc_clusters.end() );
  acc_scores.erase( remove( acc_scores.begin() , acc_scores.end() , -1.0F ) ,
                    acc_scores.end() );

}

//*****************************************************************************
// make new BasePPhoreSites at the positions of the virtual sites in sites
// so we can re-use the sites clustering functions to cluster the virtual
// sites also
void OverlayScore::extract_virtual_sites( const vector<ClusterMember> &other_sites ,
                                          vector<vector<SinglePPhoreSite *> > &virt_sites ) {

  double ddir[3];
  for( int i = 0 , is = other_sites.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "Extract virtual sites for " << other_sites[i].second->get_full_name() << endl;
#endif
    virt_sites.push_back( vector<SinglePPhoreSite *>() );
    const double *virt_site_cds = other_sites[i].second->virt_sites();
    for( int j = 0 , js = other_sites[i].second->num_virt_sites() ; j < js ; ++j ) {
      double vsite_cds[3];
      vsite_cds[0] = other_sites[i].second->coords()[0] + virt_site_cds[3*j];
      vsite_cds[1] = other_sites[i].second->coords()[1] + virt_site_cds[3*j+1];
      vsite_cds[2] = other_sites[i].second->coords()[2] + virt_site_cds[3*j+2];
      SinglePPhoreSite *vsite = new SinglePPhoreSite( vsite_cds ,
          ddir , other_sites[i].second->get_type_code() , other_sites[i].second->get_type_string() ,
          "Virtual" , false , "No Molecule" );
      virt_sites.back().push_back( vsite );
    }
  }

}

//*****************************************************************************
void OverlayScore::calc_cluster_centroid( const SiteCluster &other_sites ,
                                          double clus_centroid[3] ) {

  clus_centroid[0] = clus_centroid[1] = clus_centroid[2] = 0.0;

  for( int i = 0 , is = other_sites.second.size() ; i < is ; ++i ) {
    clus_centroid[0] += other_sites.second[i].second->coords()[0];
    clus_centroid[1] += other_sites.second[i].second->coords()[1];
    clus_centroid[2] += other_sites.second[i].second->coords()[2];
  }

  double nsites( other_sites.second.size() );
  clus_centroid[0] /= nsites;
  clus_centroid[1] /= nsites;
  clus_centroid[2] /= nsites;

}

//*****************************************************************************
float OverlayScore::mean_square_distance_to_centroid( const vector<ClusterMember> &other_sites ) {

  double centroid[3] = { 0.0 , 0.0 , 0.0 };
  BOOST_FOREACH( ClusterMember site , other_sites ) {
    const double *cds = site.second->coords();
    centroid[0] += cds[0];
    centroid[1] += cds[1];
    centroid[2] += cds[2];
  }
  double nsites = double( other_sites.size() );
  centroid[0] /= nsites;
  centroid[1] /= nsites;
  centroid[2] /= nsites;

  double dist = 0.0;
  BOOST_FOREACH( ClusterMember site , other_sites ) {
    dist += site.second->square_distance( centroid );
  }

  return float( dist / nsites );

}

//*****************************************************************************
float OverlayScore::mean_normal_cosine( const vector<ClusterMember> &other_sites ) {

  float cp = 0.0F;
  // it's not guaranteed that all sites in cluster will have a direction
  float num_sites_in_mean = 0.0F;

  for( int i = 0 , is = other_sites.size() - 1 ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "site i = " << i << " num_dirs = " << other_sites[i].second->get_num_dirs() << endl;
    other_sites[i].second->brief_report( cout );
#endif
    if( !other_sites[i].second->get_num_dirs() ) {
      continue;
    }
    const double *d1 = other_sites[i].second->direction( 0 );
    double ld1 = DACLIB::length( d1 );
    for( int j = i + 1 , js = other_sites.size() ; j < js ; ++j ) {
#ifdef NOTYET
      cout << "site j = " << j << " num_dirs = " << other_sites[j].second->get_num_dirs() << endl;
      other_sites[j].second->brief_report( cout );
#endif
      if( !other_sites[j].second->get_num_dirs() ) {
        continue;
      }
      const double *d2 = other_sites[j].second->direction( 0 );
      double ld2 = DACLIB::length( d2 );
      // the normals might be pointing in the same or opposite directions
      float tc = fabs( DACLIB::cos_angle( d1 , ld1 , d2 , ld2 ) );
      cp += tc;
      num_sites_in_mean += 1.0F;
    }
  }

  return num_sites_in_mean ? cp / num_sites_in_mean : 0.0F;

}

// ******************************************************************************
void OverlayScore::write_scores_to_stream( ostream &s , const string &sep ) const {

  s << num_sites_ << sep << rms_ << sep
    << hphobe_score_ << sep << hbond_score_ << sep << vol_score_ << sep
    << included_vol_ << sep
    << total_vol_ << sep << grid_shape_tanimoto_ << sep
    << gauss_shape_tanimoto_ << sep
    << surface_vol_ << sep << protein_clash_ << sep << mmff_nrg_ << sep
    << clique_tanimoto_ << sep << clip_score_ << sep<< num_sims_
    << sep << grid_vols_.size();
  for( int i = 0 , is = grid_vols_.size() ; i < is ; ++i ) {
    s << sep << grid_vols_[i].second;
  }

}

// ******************************************************************************
void OverlayScore::write_scores_headers_to_stream( ostream &s , const string &sep ) {

  s << "Hit_name" << sep << "Clique_size" << sep << "RMS" << sep
    << "Robin_HPhobe" << sep << "Robin_HBond" << sep << "Robin_Vol" << sep
    << "Inc.Vol" << sep << "Tot.Vol." << sep << "Grid_Shape_Tani." << sep
    << "Gauss_Shape_Tani." << sep
    << "Surf.Vol." << sep << "Prot.Clash" << sep << "MMFF_NRG" << sep
    << "Clique_Tanimoto" << sep << "Clip_Score" << sep << "Num.Sim.";
  for( int i = 0 , is = grid_vols_.size() ; i < is ; ++i ) {
    if( string( "MOE_Volume_" ) == grid_vols_[i].first.substr( 0 , 11 ) ) {
      s << sep << grid_vols_[i].first;
    } else {
      boost::filesystem::path vp( grid_vols_[i].first );
      s << sep << vp.filename();
    }
  }

}

// ******************************************************************************
// turn the scores into a set of string pairs, the first being the label,
// the second being the value.
void OverlayScore::put_scores_in_vector( vector<pair<string,string> > &scores ) const {

  scores.push_back( make_pair( string( "Grappel Num Sites" ) ,
                               lexical_cast<string>( num_sites_ ) ) );
  scores.push_back( make_pair( string( "Grappel RMS" ) ,
                               lexical_cast<string>( rms_ ) ) );
  scores.push_back( make_pair( string( "Grappel Inc. Vol." ) ,
                               lexical_cast<string>( included_vol_ ) ) );
  scores.push_back( make_pair( string( "Grappel Tot. Vol." ) ,
                               lexical_cast<string>( total_vol_ ) ) );
  if( grid_shape_tanimoto_ > -numeric_limits<float>::max() ) {
    scores.push_back( make_pair( string( "Grappel Grid Shape Tani." ) ,
                                 lexical_cast<string>( grid_shape_tanimoto_ ) ) );
  }
  if( gauss_shape_tanimoto_ > -numeric_limits<float>::max() ) {
    scores.push_back( make_pair( string( "Grappel Gauss Shape Tani." ) ,
                                 lexical_cast<string>( gauss_shape_tanimoto_ ) ) );
  }
  scores.push_back( make_pair( string( "Grappel Surf. Vol." ) ,
                               lexical_cast<string>( surface_vol_ ) ) );
  scores.push_back( make_pair( string( "Grappel Prot. Clash" ) ,
                               lexical_cast<string>( protein_clash_ ) ) );
  scores.push_back( make_pair( string( "Grappel Num. Sim." ) ,
                               lexical_cast<string>( num_sims_ ) ) );
  for( int i = 0 , is = grid_vols_.size() ; i < is ; ++i ) {
    boost::filesystem::path vp( grid_vols_[i].first );
    scores.push_back( make_pair( string( "Grappel " ) + vp.filename().string() ,
                                 lexical_cast<string>( grid_vols_[i].second ) ) );
  }

  scores.push_back( make_pair( string( "Grappel Clip Score" ) ,
			       lexical_cast<string>( clip_score_ ) ) );
  scores.push_back( make_pair( string( "Grappel MMFF" ) ,
			       lexical_cast<string>( mmff_nrg_ ) ) );
  
}

// *************************************************************************
// for packing into PVM, for example.
string OverlayScore::write_contents_to_string() {

  ostringstream oss;
  oss << "<OverlayScore>" << endl;
  oss << "<Details>" << endl;
  oss << moving_mol_name_ << endl << moving_conf_ << endl
      << fixed_mol_name_ << endl << fixed_conf_ << endl;
  oss << num_sites_ << " " << rms_
      << " " << hphobe_score_ << " " << hbond_score_ << " " << vol_score_ << " " << included_vol_
      << " " << total_vol_ << " " << grid_shape_tanimoto_ << " " << gauss_shape_tanimoto_
      << " " << surface_vol_ << " " << protein_clash_ << " " << clique_tanimoto_
      << " " << clip_score_ << " " << num_sims_ << " " << mmff_nrg_;
  oss << endl;

  oss << grid_vols_.size();
  for( int i = 0 , is = grid_vols_.size() ; i < is ; ++i ) {
    oss << " " << grid_vols_[i].first << " " << grid_vols_[i].second;
  }
  oss << endl;

  copy( site_nums_ , site_nums_ + 2 * num_sites_ ,
        ostream_iterator<int>( oss , " " ) );
  oss << endl;
  oss << "</Details>" << endl;

  if( ov_conf_ ) {
    oss << "<Molecule>" << endl;
    oemolostream oms;
    oms.SetFormat( OEFormat::SDF );
    oms.openstring();
    oms << *ov_conf_;
    oss << oms.GetString() << endl
        << "</Molecule>" << endl;
  }

  oss << "</OverlayScore>" << endl;
  
  return oss.str();

}

// ******************************************************************************
// strict weak ordering, required when sorting, needs equivalence to return false
// Effective STL, Meyers, Item 21.
bool OverlayScore::operator<( const OverlayScore &rhs ) const {

  switch( get_score_method() ) {
  case GtplDefs::GRID_SHAPE_TANI :
    return grid_shape_tani_less( rhs );
  case GtplDefs::GAUSS_SHAPE_TANI :
    return gauss_shape_tanimoto_less( rhs );
  case GtplDefs::RMS_AND_SIZE :
    return rms_and_size_less( rhs );
  case GtplDefs::SURFACE_OVLP_VOLUME :
    return surface_ovlp_volume_less( rhs );
  case GtplDefs::PROTEIN_CLASH :
    return protein_clash_less( rhs );
  case GtplDefs::MMFF_NRG :
    return mmff_nrg_less( rhs );
  case GtplDefs::ROBINS_SCORE_PARETO :
    return robins_pareto_less( rhs );
  case GtplDefs::OVERALL_SCORE_PARETO :
    return overall_pareto_less( rhs );
  }

  return false;

}

// ******************************************************************************
bool OverlayScore::operator>( const OverlayScore &rhs ) const {

  return ( rhs < *this );

}

// ******************************************************************************
bool OverlayScore::operator==( const OverlayScore &rhs ) const {

  ostringstream oss1 , oss2;
  oss1 << *this;
  oss2 << rhs;
  return( oss1.str() == oss2.str() );

}

// ******************************************************************************
bool OverlayScore::grid_shape_tani_less( const OverlayScore &rhs ) const {

  if( fabs( grid_shape_tanimoto() - rhs.grid_shape_tanimoto() ) < GtplDefs::FLT_TOL ) {
    return get_moving_mol_name() < rhs.get_moving_mol_name();
  } else {
    return( grid_shape_tanimoto() < rhs.grid_shape_tanimoto() );
  }

}

// ******************************************************************************
bool OverlayScore::gauss_shape_tanimoto_less( const OverlayScore &rhs ) const {

  if( fabs( gauss_shape_tanimoto() - rhs.gauss_shape_tanimoto() ) < GtplDefs::FLT_TOL ) {
    return get_moving_mol_name() < rhs.get_moving_mol_name();
  } else {
    return( gauss_shape_tanimoto() < rhs.gauss_shape_tanimoto() );
  }

}

// ******************************************************************************
bool OverlayScore::rms_and_size_less( const OverlayScore &rhs ) const {

  if( num_sites() == rhs.num_sites() && fabs( rms() - rhs.rms() ) < GtplDefs::FLT_TOL ) {
    return get_moving_mol_name() < rhs.get_moving_mol_name();
  } else {
    if( num_sites() == rhs.num_sites() ) {
      // bigger rms is worse
      return( rms() > rhs.rms() );
    } else {
      return( num_sites() < rhs.num_sites() );
    }
  }

}

// ******************************************************************************
bool OverlayScore::surface_ovlp_volume_less( const OverlayScore &rhs ) const {

  if( fabs( surface_volume() - rhs.surface_volume() ) < GtplDefs::FLT_TOL ) {
    return get_moving_mol_name() < rhs.get_moving_mol_name();
  } else {
    return( surface_volume() < rhs.surface_volume() );
  }

}

// ******************************************************************************
// bigger protein clash is worse
bool OverlayScore::protein_clash_less(const OverlayScore &rhs) const {

  if( fabs( protein_clash() - rhs.protein_clash() ) < GtplDefs::FLT_TOL ) {
    return get_moving_mol_name() < rhs.get_moving_mol_name();
  } else {
    return( protein_clash() > rhs.protein_clash() );
  }

}

// ******************************************************************************
bool OverlayScore::mmff_nrg_less(const OverlayScore &rhs) const {

  // bigger energy is worse
  if( fabs( get_mmff_nrg() - rhs.get_mmff_nrg() ) < GtplDefs::FLT_TOL ) {
    return get_moving_mol_name() < rhs.get_moving_mol_name();
  } else {
    return( get_mmff_nrg() > rhs.get_mmff_nrg() );
  }

}

// ******************************************************************************
bool OverlayScore::robins_pareto_less( const OverlayScore &rhs ) const {

  float lhs_score = calc_robins_pareto_score( *this );
  float rhs_score = calc_robins_pareto_score( rhs );

#ifdef NOTYET
  cout << ref_ov_->hphobe_score() << " : " << ref_ov_->hbond_score() << " : "
       << ref_ov_->vol_score() << endl;
  cout << rhs.hphobe_score() << " : " << rhs.hbond_score() << " : "
       << rhs.vol_score() << endl;
  cout << hphobe_score() << " : " << hbond_score() << " : "
       << vol_score() << endl;
  cout << "left score : " << lhs_score << " right score : " << rhs_score << endl;
#endif

  if( fabs( lhs_score - rhs_score ) < GtplDefs::FLT_TOL ) {
    return get_moving_mol_name() < rhs.get_moving_mol_name();
  } else {
    return( lhs_score < rhs_score );
  }

}

// ******************************************************************************
bool OverlayScore::overall_pareto_less( const OverlayScore &rhs ) const {

#ifdef NOTYET
  cout << "overall_pareto_less for " << get_moving_mol_name() << " vs " << get_moving_mol_name() << endl;
  cout << num_sites() << " vs " << rhs.num_sites() << " and " << rms() << " vs " << rhs.rms() << endl;
  cout << hbond_score() << " vs " << rhs.hbond_score() << endl
       << hphobe_score() << " vs " << rhs.hphobe_score() << endl
       << vol_score() << " vs " << rhs.vol_score() << endl
       << included_vol() << " vs " << rhs.included_vol() << endl
       << total_vol() << " vs " << rhs.total_vol() << endl
       << grid_shape_tanimoto() << " vs " << rhs.grid_shape_tanimoto() << endl
       << surface_volume() << " vs " << rhs.surface_volume() << endl
       << protein_clash() << " vs " << rhs.protein_clash() << endl
       << get_mmff_nrg() << " vs " << rhs.get_mmff_nrg() << endl
       << clip_score() << " vs " << rhs.clip_score() << endl
       << clique_tanimoto() << " vs " << rhs.clique_tanimoto() << endl;
#endif

  float lhs_score = calc_overall_pareto_score( *this );
  float rhs_score = calc_overall_pareto_score( rhs );

#ifdef NOTYET
  cout << "lhs_score : " << lhs_score << " and rhs_score : " << rhs_score << endl;
#endif

  if( fabs( lhs_score - rhs_score ) < GtplDefs::FLT_TOL ) {
    return get_moving_mol_name() < rhs.get_moving_mol_name();
  } else {
    return( lhs_score < rhs_score );
  }

}

// ******************************************************************************
ostream &operator<<( ostream &s , const OverlayScore &os ) {

  s << os.moving_mol_name_ << "_Conf" << os.moving_conf_ << " onto "
    << os.fixed_mol_name_ << "_Conf" << os.fixed_conf_
    << " using ";
  for( int i = 0 ; i < os.num_sites_ ; ++i ) {
    s << "(" << os.site_nums_[2*i] << "," << os.site_nums_[2*i+1] << ") ";
  }
  s << " :: " << os.rms_
    << " " << os.num_sites_ << " " << os.included_vol_ << " " << os.total_vol_
    << " " << os.grid_shape_tanimoto_
    << " " << os.gauss_shape_tanimoto_ << " " << os.surface_vol_
    << " " << os.protein_clash_;

  return s;

}

// ******************************************************************************
// check that the two molecules and the sites are the same
bool overlay_score_names_and_sites_match( const OverlayScore *a ,
                                          const OverlayScore *b ) {

  // check names match
  if( a->get_fixed_mol_name() != b->get_fixed_mol_name() ) {
    return false;
  }
  if( a->get_moving_mol_name() != b->get_moving_mol_name() ) {
    return false;
  }
  
  // check sites
  if( a->num_sites() != b->num_sites() ||
      !equal( a->sites() , a->sites() + 2 * a->num_sites() , b->sites() ) ) {
    return false;
  }

  return true;
  
}

// ********************************************************************
string remove_conf_num_from_mol_name( const string &mol_name ) {

  string ret( mol_name );
  size_t us_pos = ret.rfind( "_" );
  if( string::npos == us_pos ) {
    return ret;
  } else {
    return ret.substr( 0 , us_pos );
  }

}
