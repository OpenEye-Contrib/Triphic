//
// file PluralityHit.H
// David Cosgrove
// AstraZeneca
// 18th September 2007
//
// This class holds a hit from Plurality v3.0.

#ifndef DAC_PLURALITY_HIT
#define DAC_PLURALITY_HIT

#include "PPhoreQuery.H"
#include "VolumeGrid.H"

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

// ******************************************************************************

namespace OEChem {
  class OEMol;
  class OEMolBase;
}

class SinglePPhoreSite;

class PluralityHit {

public :

  PluralityHit() : ov_rms_( 0.0F ) , prot_overlap_( 0.0F ) ,
		   soft_exc_vol_( 0.0F ) {}
  PluralityHit( OEChem::OEMol &target_mol , int conf_num ,
                std::vector<SinglePPhoreSite *> &target_sites ,
                std::vector<int> &hit_sites ,
                boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                boost::shared_ptr<DACLIB::VolumeGrid> &soft_exc_vol_grid ,
                std::vector<BTV> &hard_exc_vols ,
                const std::vector<std::pair<std::string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                std::vector<BTC> &query_coords ,
                std::vector<BTD> &query_dists , const std::vector<float> &hit_angles ,
                const std::vector<float> &hit_torsions );
  ~PluralityHit();

  const boost::shared_ptr<OEChem::OEMolBase> get_ov_conf() const { return ov_conf_; }
  float get_ov_rms() const { return ov_rms_; }
  float get_protein_overlap() const { return prot_overlap_; }
  float get_soft_exc_vol() const { return soft_exc_vol_; }
  const std::vector<boost::shared_ptr<SinglePPhoreSite> > &get_ov_sites() const { return ov_sites_; }
  void write_to_stream( std::ostream &os , const std::string &sep );
  bool has_same_sites( PluralityHit &ph ) const;
  // add scores to the ov_conf_ as SDF tagged data
  void add_scores_to_ov_conf( const std::vector<std::string> &grid_vol_files );

  void send_via_mpi( int dest_rank );
  void receive_via_mpi( int source_rank );

private :

  boost::shared_ptr<OEChem::OEMolBase> ov_conf_;
  std::vector<boost::shared_ptr<SinglePPhoreSite> > ov_sites_;
  float ov_rms_;
  float prot_overlap_;
  float soft_exc_vol_;
  std::vector<float> hit_dists_; // just for the sites mentioned in query_dists
  std::vector<float> hit_angles_; // just for the sites mentioned in query_dists
  std::vector<float> hit_torsions_;
  std::vector<float> grid_vol_scores_; // overlap between hit and score_vol_grids
  // make the overlay conformation, adding the hit site names as tagged info
  void make_overlay_conf( OEChem::OEMol &target_mol , int conf_num );
  // overlay ov_conf_ and ov_sites_ using the transformation that moves
  // targ_cds onto query_coords.
  void overlay_hit( std::vector<BTC> &query_coords , std::vector<float> &targ_cds );

  // make the distances for sites in the query
  void fill_distances( std::vector<BTD> &query_dists );

  // calculate the overlap volume between ov_conf_ and the given grid. Creates
  // mol_grid if it doesn't already exist. Returns float::max() if it can't
  // do any better.
  float calc_overlap_vol( DACLIB::VolumeGrid &vol_grid ,
                          boost::shared_ptr<DACLIB::VolumeGrid> &mol_grid );

  // if any atom in ov_conf_ is within touching distance (distance apart less than
  // sum of 2 radii) of an excluded vol sphere, return false as the hit fails.
  bool test_hard_exc_vols( std::vector<BTV> &hard_exc_vols );

};

#endif
