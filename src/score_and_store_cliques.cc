//
// file score_and_store_cliques.cc
// David Cosgrove
// AstraZeneca
// 31st July 2007
//
// Functions to take a set of overlays, do all the necesary scoring calculations
// and put the results into a list of OverlayScore objects.

#include <list>
#include <map>
#include <string>
#include <vector>

#include <oechem.h>
#include <oeszybki.h>

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include "BasePPhoreSite.H"
#include "SinglePPhoreSite.H"
#include "OverlayScore.H"
#include "OverlayTrans.H"
#include "VolumeGrid.H"

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;
using namespace OESz;

// in eponymous file
OEMolBase *get_given_oeconf( OEMol &mol , int conf_num ,
                             bool add_conf_num_to_name );
void overlay_oemolbase( OEMolBase &mol , const OverlayTrans &ot );

// *************************************************************************
// confirm that this overlay is ok in terms of the cutoffs.
bool cutoffs_ok( const OverlayScore &ov_score ,
                 bool do_size_co , int size_cutoff ,
                 bool do_rms_co , float rms_cutoff ,
                 bool do_grid_shape_tani_co , float grid_shape_tani_cutoff ,
                 bool do_gauss_shape_tani_co , float gauss_shape_tani_cutoff ,
                 bool do_surf_ovlp_co , float surface_overlap_cutoff ,
                 bool do_inc_vol_co , float inc_vol_cutoff ,
                 bool do_prot_clash_co , float prot_clash_cutoff ,
                 bool do_mmff_co , float mmff_cutoff ) {

  if( do_size_co && ov_score.num_sites() < size_cutoff ) {
    return false;
  }

  if( do_rms_co && ov_score.rms() > rms_cutoff ) {
    return false;
  }

  if( do_grid_shape_tani_co &&
      ov_score.grid_shape_tanimoto() < grid_shape_tani_cutoff ) {
    return false;
  }

  if( do_gauss_shape_tani_co &&
      ov_score.gauss_shape_tanimoto() < gauss_shape_tani_cutoff ) {
    return false;
  }

  if( do_surf_ovlp_co && ov_score.surface_volume() < surface_overlap_cutoff ) {
    return false;
  }

  if( do_inc_vol_co && ov_score.included_vol() < inc_vol_cutoff ) {
    return false;
  }

  if( do_prot_clash_co && ov_score.protein_clash() > prot_clash_cutoff ) {
    return false;
  }

  if( do_mmff_co && ov_score.get_mmff_nrg() > mmff_cutoff ) {
    return false;
  }

  return true;

}

// *************************************************************************
// confirm that the overlay has the required points
bool required_points_satisfied( const OverlayScore &ov_score ,
                                const vector<BasePPhoreSite *> &query_sites ,
                                const vector<string> &required_points_or ,
                                const vector<string> &required_points_and ) {

#ifdef NOTYET
  cout << "OR points : ";
  copy( required_points_or.begin() , required_points_or.end() , stringOut );
  cout << endl << "AND points : ";
  copy( required_points_and.begin() , required_points_and.end() , stringOut );
  cout << endl;
#endif

  // do the OR points first, returning false if ov_score doesn't contain one of the
  // required points
  int num_points_satis = 0;
  for( int i = 0 , is = required_points_or.size() ; i < is ; ++i ) {
    for( int j = 0 , js = ov_score.num_sites() ; j < js ; ++j ) {
      if( query_sites[ov_score.sites()[2*j]]->get_type_string() ==
          required_points_or[i] ) {
        num_points_satis = 1;
        break;
      }
    }
    if( num_points_satis )
      break; // we only need one of them
  }
  if( !required_points_or.empty() && !num_points_satis )
    return false;

  // now do the ANDs. If a point was in the OR list, we need two different sites of
  // that type.
  map<string,int> num_needed;
  for( int i = 0 , is = required_points_or.size() ; i < is ; ++i )
    num_needed.insert( make_pair( required_points_or[i] , 1 ) );
  map<string,int>::iterator p;
  for( int i = 0 , is = required_points_and.size() ; i < is ; ++i ) {
    p = num_needed.find( required_points_and[i] );
    if( p == num_needed.end() )
      num_needed.insert( make_pair( required_points_and[i] , 1 ) );
    else
      p->second++;
  }

  for( int i = 0 , is = required_points_and.size() ; i < is ; ++i ) {
    num_points_satis = 0;
    for( int j = 0 , js = ov_score.num_sites() ; j < js ; ++j ) {
      if( query_sites[ov_score.sites()[2*j]]->get_type_string() ==
          required_points_and[i] ) {
        ++num_points_satis;
      }
    }
    p = num_needed.find( required_points_and[i] );
    if( num_points_satis < p->second )
      return false;
  }

  return true;

}

// *************************************************************************
// confirm that the overlay has the required sites
bool required_sites_satisfied( const OverlayScore &ov_score ,
                               const vector<BasePPhoreSite *> &query_sites ,
                               const vector<SinglePPhoreSite *> &target_sites ,
                               const vector<string> &required_sites_or ,
                               const vector<string> &required_sites_and ) {

#ifdef NOTYET
  cout << "OR sites : ";
  copy( required_sites_or.begin() , required_sites_or.end() , stringOut );
  cout << endl << "AND sites : ";
  copy( required_sites_and.begin() , required_sites_and.end() , stringOut );
  cout << endl;
#endif

  // count the number of times the ORs and ANDs are satisfied
  unsigned int num_or_sites_satis = 0 , num_and_sites_satis = 0;
  for( int j = 0 , js = ov_score.num_sites() ; j < js ; ++j ) {
    string query_site_name =
        query_sites[ov_score.sites()[2*j]]->get_full_name();
    string target_site_name =
        target_sites[ov_score.sites()[2*j+1]]->get_full_name();
    for( int i = 0 , is = required_sites_or.size() ; i < is ; ++i ) {
      if( query_site_name == required_sites_or[i] ||
          target_site_name == required_sites_or[i]) {
        num_or_sites_satis = 1;
        break;
      }
    }
    for( int i = 0 , is = required_sites_and.size() ; i < is ; ++i ) {
      if( query_site_name == required_sites_and[i] ||
          target_site_name == required_sites_and[i]) {
        ++num_and_sites_satis;
      }
    }
  }
  if( !required_sites_or.empty() && !num_or_sites_satis )
    return false;
  if( num_and_sites_satis < required_sites_and.size() )
    return false;
  return true;

}

// ****************************************************************************
void add_query_clique_site_info( const vector<BasePPhoreSite *> &query_sites ,
                                 const vector<int> &clique , OEMolBase &mol ) {

  string site_names;
  for( unsigned int i = 0 , is = clique.size() ; i < is ; i += 2 )
    site_names += query_sites[clique[i]]->get_full_name() + " ";

  // chop the last space off
  site_names = site_names.substr( 0 , site_names.length() - 1 );

  OESetSDData( mol , "Query_PPhore_Site_Names" , site_names );

}

// ****************************************************************************
void add_target_clique_site_info( const vector<SinglePPhoreSite *> &target_sites ,
                                  const vector<int> &clique , OEMolBase &mol ) {

  string site_names;
  for( unsigned int i = 1 , is = clique.size() ; i < is ; i += 2 )
    site_names += target_sites[clique[i]]->get_full_name() + " ";

  // chop the last space off
  site_names = site_names.substr( 0 , site_names.length() - 1 );

  OESetSDData( mol , "Target_PPhore_Site_Names" , site_names );

}

// ****************************************************************************
void add_clique_site_info( const vector<BasePPhoreSite *> &query_sites ,
                           const vector<SinglePPhoreSite *> &target_sites ,
                           const vector<int> &clique , OEMolBase &mol ) {

  add_query_clique_site_info( query_sites , clique , mol );
  add_target_clique_site_info( target_sites , clique , mol );

}

// **************************************************************************
// non-member function that takes a multi-molecule query and corresponding
// SinglePPhoreSites and overlays them with the given clique, returning a new
// copy of target conf and sites, transformed to the overlay.
void overlay_mols_and_sites( OEMol *target_mol ,
                             const vector<BasePPhoreSite *> &query_sites ,
                             const vector<SinglePPhoreSite *> &target_sites ,
                             OverlayScore *ov_score , bool use_ring_norms ,
                             bool use_h_vectors , bool use_lps ,
                             OEMolBase *&target_conf ) {

  int target_conf_num = ov_score->get_moving_conf();
  target_conf = get_given_oeconf( *target_mol , target_conf_num , false );
  if( !target_conf ) {
    return;
  }

  vector<SinglePPhoreSite *> ov_target_sites;
  for( int j = 0 , js = target_sites.size() ; j < js ; ++j ) {
    ov_target_sites.push_back( new SinglePPhoreSite( *target_sites[j] ) );
  }

  OverlayTrans ov_trans;
  vector<int> clique( ov_score->sites() ,
                      ov_score->sites() + ov_score->num_sites() * 2 );
  calc_overlay_trans( query_sites , ov_target_sites , clique , ov_trans ,
                      false , false , false );

  if( use_ring_norms || use_h_vectors || use_lps ) {
    calc_overlay_trans( query_sites , ov_target_sites , clique , ov_trans ,
                        use_ring_norms , use_h_vectors , use_lps );
  }

  overlay_oemolbase( *target_conf , ov_trans );
  overlay( ov_target_sites , ov_trans );

  for( int i = 0 , is = clique.size() ; i < is ; i += 2 ) {
    if( ov_target_sites[clique[i+1]]->get_twiddlable() &&
        use_h_vectors ) {
      ov_target_sites[clique[i+1]]->twiddle( *query_sites[clique[i]] ,
          GtplDefs::H_VECTOR );
    }
  }
  add_clique_site_info( query_sites , ov_target_sites , clique , *target_conf );

  for( int i = 0 , is = ov_target_sites.size() ; i < is ; ++i ) {
    delete ov_target_sites[i];
  }

}

// **************************************************************************
// free function that takes SinglePPhoreSites and overlays them using the
// OverlayTrans in the OverlayScore, returning copies of the sites, transformed
// to the overlay.
void overlay_sites( const vector<BasePPhoreSite *> &query_sites ,
                    const vector<SinglePPhoreSite *> &target_sites ,
                    OverlayScore *ov_score , bool use_h_vectors ,
                    vector<SinglePPhoreSite *> &ov_target_sites ) {

  ov_target_sites.clear();

  for( int j = 0 , js = target_sites.size() ; j < js ; ++j ) {
    ov_target_sites.push_back( new SinglePPhoreSite( *target_sites[j] ) );
  }

  overlay( ov_target_sites , *ov_score->ov_trans() );

  const int *sites = ov_score->sites();
  for( int i = 0 , is = 2 * ov_score->num_sites() ; i < is ; i += 2 ) {
    if( ov_target_sites[sites[i+1]]->get_twiddlable() && use_h_vectors )
      ov_target_sites[sites[i+1]]->twiddle( *query_sites[sites[i]] ,
          GtplDefs::H_VECTOR );
  }

}

// **************************************************************************
// free function that takes a multi-molecule query overlays it using the
// OverlayTrans in the OverlayScore, returning a new copy of target conf
// to the overlay.
void overlay_mol( OEMol *target_mol ,
                  const vector<BasePPhoreSite *> &query_sites ,
                  const vector<SinglePPhoreSite *> &ov_target_sites ,
                  OverlayScore *ov_score , OEMolBase *&target_conf ) {

  int target_conf_num = ov_score->get_moving_conf();

  target_conf = get_given_oeconf( *target_mol , target_conf_num , false );
  if( !target_conf ) {
    return;
  }
  overlay_oemolbase( *target_conf , *ov_score->ov_trans() );

  vector<int> clique( ov_score->sites() ,
                      ov_score->sites() + ov_score->num_sites() * 2 );
  add_clique_site_info( query_sites , ov_target_sites , clique , *target_conf );

}

// *************************************************************************
void increment_score_ranks( vector<pair<OverlayScore *,float> > &scores ,
                            map<OverlayScore *,int> &rank_scores ,
                            bool descending_sort ) {

  if( descending_sort ) {
    sort( scores.begin() , scores.end() ,
          bind( greater<float>() ,
                bind( &pair<OverlayScore *,float>::second , _1 ) ,
                bind( &pair<OverlayScore *,float>::second , _2 ) ) );
  } else {
    sort( scores.begin() , scores.end() ,
          bind( less<float>() ,
                bind( &pair<OverlayScore *,float>::second , _1 ) ,
                bind( &pair<OverlayScore *,float>::second , _2 ) ) );
  }

  float last_score = scores.front().second;
  int curr_rank = 1;

  typedef pair<OverlayScore *,float> OSF;
  BOOST_FOREACH( OSF s , scores ) {
    map<OverlayScore *,int>::iterator p = rank_scores.find( s.first );
    if( s.second != last_score ) {
      ++curr_rank;
    }
    p->second += curr_rank;
    last_score = s.second;
  }

}

// *************************************************************************
void add_overlay_score_to_normal_list( OverlayScore *new_score ,
                                       list<OverlayScore *> &clique_overlays ) {


  // put the new score into clique_overlays
  list<OverlayScore *>::iterator where =
      lower_bound( clique_overlays.begin() , clique_overlays.end() ,
                   new_score , OverlayScoreMore() );
  clique_overlays.insert( where , new_score );

#ifdef NOTYET
  if( !clique_overlays.empty() ) {
    cout << "front score : " << clique_overlays.front()->get_moving_mol_name() << " : " << clique_overlays.front()->get_moving_conf() << " : ";
    clique_overlays.front()->write_scores_to_stream( cout , " " );
    cout << endl;
    cout << "back score : " << clique_overlays.back()->get_moving_mol_name() << " : " << clique_overlays.back()->get_moving_conf() << " : ";
    clique_overlays.back()->write_scores_to_stream( cout , " " );
    cout << endl;
    cout << "insertion position : " << std::distance( clique_overlays.begin() , where ) << endl;
  }
#endif

}

// *************************************************************************
void add_overlay_score_to_list( OverlayScore *new_score ,
                                list<OverlayScore *> &clique_overlays ) {

#ifdef NOTYET
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  cout << "This hit : " << new_score->get_moving_mol_name() << " : " << new_score->get_moving_conf() << " : ";
  // new_score->write_scores_to_stream( cout , " " );
  cout << new_score->hphobe_score() << " , " << new_score->hbond_score() << " , " << new_score->vol_score();
  cout << endl;
  cout << "Current list" << endl;
  for( list<OverlayScore *>::iterator p = clique_overlays.begin() ; p != clique_overlays.end() ; ++p ) {
    cout << std::distance( clique_overlays.begin() , p ) << " : "
         << (*p)->get_moving_mol_name() << " : " << (*p)->get_moving_conf() << " : ";
    //   (*p)->write_scores_to_stream( cout , " " );
    cout << (*p)->hphobe_score() << " , " << (*p)->hbond_score() << " , " << (*p)->vol_score();
    cout << endl;
  }
#endif

  if( clique_overlays.empty() ) {
#ifdef NOTYET
    cout << "first one" << endl;
#endif
    clique_overlays.push_back( new_score );
    return;
  }

  // if this overlay already exists, add as a similar.
  list<OverlayScore *>::iterator where;
  for( where = clique_overlays.begin() ; where != clique_overlays.end() ; ++where ) {
    if( overlay_score_names_and_sites_match( *where , new_score ) ) {
      break;
    }
  }

  if( where != clique_overlays.end() ) {
#ifdef NOTYET
    cout << "Found a similar" << endl;
    cout << "where : ";
    (*where)->write_scores_to_stream( cout , " " );
    cout << endl;
    cout << "new_score : ";
    new_score->write_scores_to_stream( cout , " " );
    cout << endl;
#endif
    if( *(*where) == *new_score ) {
      // if they're the same, they're not similar. This normally arises if it's
      // a restart job and it's re-doing a bit of the search.
      delete new_score;
      return;
    }
    // if (*where) is better than new_score, add new_score to (*where), otherwise
    // need to add (*where) and all its similars to new_score
    if( *(*where) > *new_score ) {
      (*where)->add_similar();
      delete new_score; // don't need it anymore
      return;
    } else {
      // add (*where) to new_score and remove (*where) from clique_overlays
      // new_score needs to be added to clique_overlays in the appropriate
      // place, so don't return at this point.
      new_score->add_similars( *where );
      delete *where; // don't need it any more
      clique_overlays.erase( where );
    }
  }

  // pareto ranking methods need to be dealt with differently from those that
  // use a straightforward score
  add_overlay_score_to_normal_list( new_score , clique_overlays );

#ifdef NOTYET
  cout << "New list" << endl;
  for( list<OverlayScore *>::iterator p = clique_overlays.begin() ; p != clique_overlays.end() ; ++p ) {
    cout << std::distance( clique_overlays.begin() , p ) << " : "
         << (*p)->get_moving_mol_name() << " : " << (*p)->get_moving_conf() << " : ";
    cout << (*p)->hphobe_score() << " , " << (*p)->hbond_score() << " , " << (*p)->vol_score();
    cout << endl;
  }
#endif

}

// ****************************************************************************
// take the cliques for the given conformations of the given molecules, score
// them and keep the top so many in clique_overlays.
// If max_num_hits == 0 there's no maximum.
// In a slightly confusing manner, the target is being overlaid onto the query,
// i.e. query_sites are fixed, target_sites are moving.
// I'm not sure how that came about.
void score_and_store_cliques( const string &query_name , int query_conf_num ,
                              OEMolBase *query_conf ,
                              vector<BasePPhoreSite *> &query_sites ,
                              vector<SinglePPhoreSite *> &query_score_sites ,
                              vector<SinglePPhoreSite *> &target_sites ,
                              OEMol *target_mol , int target_conf_num ,
                              boost::shared_ptr<DACLIB::VolumeGrid> &query_solid_grid ,
                              boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                              const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                              const vector<vector<int> > &cliques ,
                              GtplDefs::DIRS_USAGE ring_norm_usage ,
                              GtplDefs::DIRS_USAGE h_vector_usage ,
                              GtplDefs::DIRS_USAGE lp_usage ,
                              float ring_norm_tol , float h_vector_tol ,
                              float lp_tol , bool do_size_co , int size_cutoff ,
                              bool do_rms_co , float rms_cutoff ,
                              bool do_grid_shape_tani_co ,
                              float grid_shape_tani_cutoff ,
                              bool do_gauss_shape_tani_co ,
                              float gauss_shape_tani_cutoff ,
                              bool do_surf_ovlp_co , float surface_overlap_cutoff ,
                              bool do_inc_vol_co , float inc_vol_cutoff ,
                              bool do_prot_clash_co , float prot_clash_cutoff ,
                              bool do_mmff_co , float mmff_cutoff ,
                              const vector<string> &required_sites_or ,
                              const vector<string> &required_sites_and ,
                              const vector<string> &required_points_or ,
                              const vector<string> &required_points_and ,
                              boost::shared_ptr<OESzybki> &szybki ,
                              list<OverlayScore *> &clique_overlays ) {

  string target_title = target_mol->GetTitle();

  vector<boost::shared_ptr<SinglePPhoreSite> > ov_target_sites;

  bool use_ring_norms_in_overlay = ring_norm_usage == GtplDefs::ALIGN;
  bool check_ring_norms = ( ring_norm_usage == GtplDefs::SCREEN ||
                            ring_norm_usage == GtplDefs::ALIGN );
  float rnt = cos( M_PI * ring_norm_tol / 180.0 );

  bool use_h_vectors_in_overlay = h_vector_usage == GtplDefs::ALIGN;
  bool check_h_vectors = ( h_vector_usage == GtplDefs::SCREEN ||
                           h_vector_usage == GtplDefs::ALIGN );
  float hvt = cos( M_PI * h_vector_tol / 180.0 );

  bool use_lps_in_overlay = lp_usage == GtplDefs::ALIGN;
  bool check_lps = ( lp_usage == GtplDefs::SCREEN ||
                     lp_usage == GtplDefs::ALIGN );
  float lpt = cos( M_PI * lp_tol / 180.0 );

#ifdef NOTYET
  cout << "Scoring " << cliques.size() << " cliques." << endl;
#endif

  for( int i = 0 , is = cliques.size() ; i < is ; ++i ) {

#ifdef NOTYET
    cout << "Clique " << i << " : " << cliques[i].size() / 2 << " :: ";
    for( int j = 0 , js = cliques[i].size() ; j < js ; j += 2 ) {
      cout << "(" << cliques[i][j] << "," << cliques[i][j+1] << ") ";
    }
    cout << "   Moving conf num = " << target_conf_num << endl;
#endif
    OverlayScore *new_score = 0;
    try {
      new_score = new OverlayScore( query_name , target_title , query_conf_num ,
                                    target_conf_num , cliques[i] , query_sites ,
                                    query_score_sites , target_sites ,
                                    *query_conf , *target_mol ,
                                    query_solid_grid , protein_grid , score_vol_grids ,
                                    szybki , use_ring_norms_in_overlay ,
                                    use_h_vectors_in_overlay , use_lps_in_overlay ,
                                    true /* do the overlay */ );
    } catch( OverlayScoreError &e ) {
      // it failed, most likely because there were less than 3 sites with
      // unique positions in the clique. Sometimes you get 2 sites based on the
      // same atom (e.g. donor and acceptor on a hydroxyl)
      continue;
    }
    ov_target_sites = new_score->get_ov_sites();

    // make sure any vectors are lined up, including twiddling/flipping
    // X-H or LP vectors if possible. This is in SinglePPhoreSite.cc
    if( confirm_site_vectors( cliques[i] , query_sites , ov_target_sites ,
                              check_ring_norms , rnt , check_h_vectors , hvt ,
                              check_lps , lpt ) ) {
      if( new_score->total_vol() < 0.0F ||
          do_grid_shape_tani_co || do_inc_vol_co ||
          do_surf_ovlp_co || do_prot_clash_co ) {
        new_score->calc_volume_scores( new_score->get_ov_conf() , query_solid_grid ,
                                       protein_grid );
      }
      if( new_score->gauss_shape_tanimoto() < 0.0F && do_gauss_shape_tani_co ) {
        new_score->calc_gauss_shape_tanimoto( *query_conf , *new_score->get_ov_conf() );
      }

      // if we're optimising, it'll already be done, but single point calculation
      // delayed until necessary
      if( szybki && OERunType::SinglePoint == szybki->GetRunType() &&
          numeric_limits<float>::max() != mmff_cutoff ) {
        OESzybkiResults sr;
        float mmff_nrg;
        if( !(*szybki)( *new_score->get_ov_conf() , sr ) ) {
          cerr << "Optimisation failed for " << new_score->get_ov_conf()->GetTitle() << endl;
          mmff_nrg = numeric_limits<float>::max();
        } else {
          mmff_nrg = sr.GetTotalEnergy();
        }
        new_score->set_mmff_nrg( mmff_nrg );
      }

      // we're now in a position to verify if cutoffs are satisfied.
      if( cutoffs_ok( *new_score , do_size_co , size_cutoff ,
                      do_rms_co , rms_cutoff , do_grid_shape_tani_co ,
                      grid_shape_tani_cutoff , do_gauss_shape_tani_co ,
                      gauss_shape_tani_cutoff , do_surf_ovlp_co ,
                      surface_overlap_cutoff , do_inc_vol_co , inc_vol_cutoff ,
                      do_prot_clash_co , prot_clash_cutoff ,
                      do_mmff_co , mmff_cutoff ) &&
          required_points_satisfied( *new_score , query_sites ,
                                     required_points_or , required_points_and ) &&
          required_sites_satisfied( *new_score , query_sites , target_sites ,
                                    required_sites_or , required_sites_and ) ) {
        add_overlay_score_to_list( new_score , clique_overlays );
      } else {
        delete new_score;
      }
    } else {
      delete new_score; // didn't make it
    }

  }

}

// ****************************************************************************
// take the cliques for the given conformations of the given molecules, score
// them and keep the top so many in clique_overlays
// This one, as you might surmise, doesn't move the molecules, just scores them.
void score_and_store_cliques_no_overlay( const string &query_name ,
                                         OEMolBase *query_conf ,
                                         vector<BasePPhoreSite *> &query_sites ,
                                         vector<SinglePPhoreSite *> &query_score_sites ,
                                         vector<SinglePPhoreSite *> &target_sites ,
                                         OEMol *target_mol , int target_conf_num ,
                                         boost::shared_ptr<DACLIB::VolumeGrid> &query_solid_grid ,
                                         boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                                         const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                                         const vector<vector<int> > &cliques ,
                                         GtplDefs::DIRS_USAGE ring_norm_usage ,
                                         GtplDefs::DIRS_USAGE h_vector_usage ,
                                         GtplDefs::DIRS_USAGE lp_usage ,
                                         float ring_norm_tol , float h_vector_tol ,
                                         float lp_tol , int size_cutoff ,
                                         float rms_cutoff ,
                                         bool do_grid_shape_tani_co ,
                                         float grid_shape_tani_cutoff ,
                                         bool do_gauss_shape_tani_co ,
                                         float gauss_shape_tani_cutoff ,
                                         bool do_surf_ovlp_co , float surface_overlap_cutoff ,
                                         bool do_inc_vol_co , float inc_vol_cutoff ,
                                         bool do_prot_clash_co , float prot_clash_cutoff ,
                                         const vector<string> &required_sites_or ,
                                         const vector<string> &required_sites_and ,
                                         const vector<string> &required_points_or ,
                                         const vector<string> &required_points_and ,
                                         list<OverlayScore *> &clique_overlays ) {

  string target_title = target_mol->GetTitle();

  bool use_ring_norms_in_overlay = ring_norm_usage == GtplDefs::ALIGN;
  bool check_ring_norms = ( ring_norm_usage == GtplDefs::SCREEN ||
                            ring_norm_usage == GtplDefs::ALIGN );
  float rnt = cos( M_PI * ring_norm_tol / 180.0 );

  bool use_h_vectors_in_overlay = h_vector_usage == GtplDefs::ALIGN;
  bool check_h_vectors = ( h_vector_usage == GtplDefs::SCREEN ||
                           h_vector_usage == GtplDefs::ALIGN );
  float hvt = cos( M_PI * h_vector_tol / 180.0 );

  bool use_lps_in_overlay = lp_usage == GtplDefs::ALIGN;
  bool check_lps = ( lp_usage == GtplDefs::SCREEN ||
                     lp_usage == GtplDefs::ALIGN );
  float lpt = cos( M_PI * lp_tol / 180.0 );

  boost::shared_ptr<OESzybki> szybki; // not used

  for( int i = 0 , is = cliques.size() ; i < is ; ++i ) {

    OverlayScore *new_score = 0;
    try {
      new_score = new OverlayScore( query_name , target_title , 0 , target_conf_num ,
                                    cliques[i] , query_sites , query_score_sites ,
                                    target_sites , *query_conf , *target_mol ,
                                    query_solid_grid , protein_grid , score_vol_grids ,
                                    szybki , use_ring_norms_in_overlay ,
                                    use_h_vectors_in_overlay , use_lps_in_overlay ,
                                    false /* don't do overlay */);
    } catch( OverlayScoreError &e ) {
      // it failed, most likely because there were less than 3 sites with
      // unique positions in the clique. Sometimes you get 2 sites based on the
      // same atom (e.g. donor and acceptor on a hydroxyl)
      continue;
    }

    // make sure any vectors are lined up, including twiddling/flipping
    // X-H or LP vectors if possible. This is in SinglePPhoreSite.cc
    vector<boost::shared_ptr<SinglePPhoreSite> > ov_target_sites = new_score->get_ov_sites();
    if( confirm_site_vectors( cliques[i] , query_sites , ov_target_sites ,
                              check_ring_norms , rnt , check_h_vectors , hvt ,
                              check_lps , lpt ) ) {
      if( new_score->total_vol() < 0.0F ||
          do_grid_shape_tani_co || do_inc_vol_co ||
          do_surf_ovlp_co || do_prot_clash_co ) {
        // do the volume scores
        new_score->calc_volume_scores( new_score->get_ov_conf() , query_solid_grid , protein_grid );
      }
      if( new_score->gauss_shape_tanimoto() < 0.0F && do_gauss_shape_tani_co ) {
        new_score->calc_gauss_shape_tanimoto( *query_conf , *new_score->get_ov_conf() );
      }

      // we're now in a position to verify if cutoffs are satisfied.
      if( cutoffs_ok( *new_score , true , size_cutoff ,
                      true , rms_cutoff , do_grid_shape_tani_co ,
                      grid_shape_tani_cutoff , do_gauss_shape_tani_co ,
                      gauss_shape_tani_cutoff , do_surf_ovlp_co ,
                      surface_overlap_cutoff , do_inc_vol_co , inc_vol_cutoff ,
                      do_prot_clash_co , prot_clash_cutoff ,
                      false , -1.0 ) &&
          required_points_satisfied( *new_score , query_sites ,
                                     required_points_or , required_points_and ) &&
          required_sites_satisfied( *new_score , query_sites , target_sites ,
                                    required_sites_or , required_sites_and ) ) {
        add_overlay_score_to_list( new_score , clique_overlays );
      } else {
        delete new_score;
      }
    } else {
      delete new_score; // didn't make it
    }

  }

}

