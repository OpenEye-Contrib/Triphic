//
// file OverlayScore.H
// David Cosgrove
// AstraZeneca
// 21st June 2006
//
// This files is the interface file for the class OverlayScore, which holds the
// details of a clique based overlay between 2 conformations of 2 molecules
// and the scores for the overlay.

#ifndef DAC_OVERLAY_SCORE__
#define DAC_OVERLAY_SCORE__

#include <functional>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "gtpl_defs.H"
#include "VolumeGrid.H"

// ******************************************************************************

class BasePPhoreSite;
class OverlayTrans;
class SinglePPhoreSite;

namespace OEChem {
class OEMolBase;
}

namespace OESz {
class OESzybki;
}

class OverlayScoreError {

public :
  explicit OverlayScoreError( const std::string &msg ) : msg_( msg ) {}
  const char *what() const { return msg_.c_str(); }
private :
  std::string msg_;

};

class OverlayScore {

public :

  // These typedefs are used for calculating Robin's hphobe, donor and acceptor
  // scores.
  // ClusterMember is a molecule number and a site
  typedef std::pair<int,SinglePPhoreSite *> ClusterMember;
  // cluster is a set of ClusterMembers and a mean square distance between members
  // the mean square distance will only be calculated if needed
  typedef std::pair<float,std::vector<ClusterMember> > SiteCluster;

protected :

  // the names of the two molecules
  std::string fixed_mol_name_ , moving_mol_name_;
  int    fixed_conf_ , moving_conf_; // conformations of those mols.

  float        rms_;
  int          num_sites_; // in the overlay
  int          site_nums_[100]; /* paired, as they came from the clique, so
           clique can't be bigger than 50 pairs */

  // Robin's overlay scores from JCAMD 26, 451-472(2012).
  // vol_score_ should be the same as total_vol_ but keeping them consistent
  // for now.
  float        hphobe_score_ , hbond_score_ , vol_score_;

  // volume-based scores
  float        included_vol_; // intersection volume of all mols
  float        total_vol_; // total occupied volume
  float        grid_shape_tanimoto_; /* calculated between first n-1 mols and nth,
					using the VolumeGrid method */
  float        gauss_shape_tanimoto_; /* calculated using ShapeTanimoto object,
					 only used in Triphic but more convenient
					 here. bad design, I guess, but this is
					 a retro-fit. */
  float        surface_vol_; // included volumes of all molecules
  float        protein_clash_; /* included volume with another molecule, presumed
				  to be a protein active site. -FLOAT_MAX if
				  not used. */
  float        mmff_nrg_; // between protein and ligand
  int          num_sims_; /* number of overlays with the same sites
			     as this one, but with poorer scores */
  float        clip_score_ , clique_tanimoto_;


  // intersection volume with 1 or more named arbitrary VolumeGrids.
  std::vector<std::pair<std::string,float> > grid_vols_;

  // score method for comparison between 2 scores
  static GtplDefs::SCORE_METHOD score_method_;

  boost::scoped_ptr<OverlayTrans> ov_trans_;
  // copies of the moving sites and conformation, overlaid using the clique
  boost::shared_ptr<OEChem::OEMolBase> ov_conf_;
  std::vector<boost::shared_ptr<SinglePPhoreSite> > ov_sites_;

  void copy_data( const OverlayScore &ov );

  // volume scores off a grid
  void calc_volume_scores( boost::scoped_ptr<DACLIB::VolumeGrid> &mol_grid ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &target_solid_grid ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid );

  void optimise_overlay( boost::shared_ptr<OESz::OESzybki> &szybki );

public :

  OverlayScore();
  OverlayScore( const std::string &fix_mol_name , const std::string &mov_mol_name ,
                int fix_conf_num , int mov_conf_num , const std::vector<int> &clique ,
                const std::vector<BasePPhoreSite *> &fixed_sites ,
                const std::vector<std::vector<SinglePPhoreSite *> > &fixed_score_sites ,
                const std::vector<SinglePPhoreSite *> &moving_sites ,
                OEChem::OEMolBase &fixed_conf , OEChem::OEMol &moving_mol ,
                boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid ,
                boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                const std::vector<std::pair<std::string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                boost::shared_ptr<OESz::OESzybki> &szybki ,
                bool use_dirs_in_overlay , bool use_h_vectors_in_overlay ,
                bool use_lps_in_overlay , bool do_overlay = true );

  // for the comparisons
  OverlayScore( const std::string &fix_mol_name , const std::string &mov_mol_name );

  // this one builds itself from a string, and is used, at least initially,
  // in parallel triphic to create an object from data sent between
  // master and slave.
  OverlayScore( const std::string &str );

  OverlayScore( const OverlayScore &ovs );

  virtual ~OverlayScore();

  // make a copy of the this OverlayScore, apart from the similars and the
  // OverlayTrans
  OverlayScore *make_copy_no_sims();

  // increment the number of similar overlays.
  void add_similar( OverlayScore *sim );
  // copy the similars of sims to this, remove them from sims, and add sims to this
  void add_similars( OverlayScore *sims );
  // clear the similars
  void clear_similars();

  std::string get_fixed_mol_name() const { return fixed_mol_name_; }
  std::string get_moving_mol_name() const { return moving_mol_name_; }
  int get_fixed_conf() const { return fixed_conf_; }
  int get_moving_conf() const { return moving_conf_; }
  void set_ov_conf( OEChem::OEMolBase *new_ov_conf );
  OEChem::OEMolBase *get_ov_conf();
  std::vector<boost::shared_ptr<SinglePPhoreSite> > get_ov_sites() const { return ov_sites_; }
  void set_ov_sites( std::vector<SinglePPhoreSite *> &new_ov_sites );

  float rms() const { return rms_; }
  void set_rms( float new_rms ) { rms_ = new_rms; }
  void set_mmff_nrg( float new_val ) { mmff_nrg_ = new_val; }
  float get_mmff_nrg() const { return mmff_nrg_; }

  int num_sites() const { return num_sites_; }
  const int *sites() const { return site_nums_; }
  int num_sim() const { return num_sims_; }
  float included_vol() const { return included_vol_; }
  float total_vol() const { return total_vol_; }
  float grid_shape_tanimoto() const { return grid_shape_tanimoto_; }
  float gauss_shape_tanimoto() const { return gauss_shape_tanimoto_; }
  float surface_volume() const { return surface_vol_; }
  float protein_clash() const { return protein_clash_; }
  float clip_score() const { return clip_score_; }
  float clique_tanimoto() const { return clique_tanimoto_; }
  
  float hphobe_score() const { return hphobe_score_; }
  float hbond_score() const { return hbond_score_; }
  float vol_score() const { return vol_score_; }

  OverlayTrans *ov_trans() const { return ov_trans_.get(); }
  int num_grid_vols() const { return grid_vols_.size(); }
  std::vector<std::pair<std::string,float> > grid_vols() const {
    return grid_vols_;
  }
  void add_grid_vol( std::string &nam , float vol ) {
    grid_vols_.push_back( std::make_pair( nam , vol ) );
  }

  // calculate the volume based scores
  void calc_gauss_shape_tanimoto( OEChem::OEMolBase &mol1 ,
                                  OEChem::OEMolBase &mol2 );
  void calc_volume_scores( OEChem::OEMolBase *mol ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &target_solid_grid ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid );
  void calc_volume_scores( OEChem::OEMolBase *mol ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &target_solid_grid ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                           const std::vector<std::pair<std::string,DACLIB::VolumeGrid *> > &vol_grids );
  void calc_protein_clash( boost::shared_ptr<DACLIB::VolumeGrid> &mol_grid ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &prot_grid );

  void overlay_mol_and_sites( const std::vector<BasePPhoreSite *> &fixed_sites ,
                              const std::vector<SinglePPhoreSite *> &moving_sites ,
                              OEChem::OEMol &moving_mol ,
                              boost::shared_ptr<OESz::OESzybki> &szybki,
                              bool use_ring_norms_in_overlay ,
                              bool use_h_vectors_in_overlay , bool use_lps_in_overlay ,
                              bool do_overlay );

  // overlay moving sites using ov_trans_ to make ov_sites_. If use_h_vector or
  // use_lp_vector are true, the relevant directions are altered to have them
  // line up with the query
  void overlay_sites( const std::vector<BasePPhoreSite *> &fixed_sites ,
                      const std::vector<SinglePPhoreSite *> &moving_sites ,
                      bool use_h_vectors , bool use_lp_vectors );
  // free function that takes a multi-molecule query overlays it using the
  // OverlayTrans in the OverlayScore, returning a new copy of target conf
  // to the overlay.
  void overlay_moving_conf( OEChem::OEMol &target_mol ,
                            const std::vector<BasePPhoreSite *> &fixed_sites ,
                            const std::vector<SinglePPhoreSite *> &moving_sites);
  void add_clique_fixed_site_info( const std::vector<BasePPhoreSite *> &fixed_sites ,
                                   OEChem::OEMolBase &mol );
  void add_clique_moving_site_info( const std::vector<SinglePPhoreSite *> &moving_sites ,
                                    OEChem::OEMolBase &mol );
  // clique tanimoto is just the proportion of search sites that appear in the
  // clique. Not sure anyone ever uses it, but triphic has always calculated it.
  void calc_clique_tanimoto( int num_query_sites , int num_target_sites );
  // calculate the clip score. Clip is a program similar to Triphic that Willet
  // et al have published (JCICS,43,443-448 (2003)). The original Triphic
  // pre-dated the publication.
  // The clip score in this case is the clique tanimoto
  // score times ( 1.0 - RMSD( matching distances ) / ( sum of max poss distance
  // error )
  void calc_clip_score( std::vector<BasePPhoreSite *> &query_sites ,
                        float dist_tol );
  
  // calculate the hphobe, donor and acceptor scores. The target volume grid is
  // needed for calculating the occlusion of virtual donor and acceptor sites.
  void calc_robins_scores( const std::vector<std::vector<SinglePPhoreSite *> > &fixed_score_sites ,
                           const std::vector<SinglePPhoreSite *> &moving_sites ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid );
  void calc_directional_hphobe_score( const std::vector<std::vector<SinglePPhoreSite *> > &sites );
  void calc_hbond_score( const std::vector<std::vector<SinglePPhoreSite *> > &sites ,
                         boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid );
  // cluster sites of given type not caring about whether their possession of a
  // direction. Only 1 site in each cluster from each molecule
  void cluster_sites( const std::string &site_type , float thresh ,
                      const std::vector<std::vector<SinglePPhoreSite *> > &sites ,
                      std::vector<SiteCluster> &clusters );
  // make near-neighbour lists of sites with type string site_type
  void make_site_nnls( const std::string &site_type , float thresh ,
                       const std::vector<std::vector<SinglePPhoreSite *> > &sites ,
                       std::vector<SiteCluster> &nnls );
  // calculate the sum of square distances of the cluster/nnl
  void calculate_nnl_dists( SiteCluster &nnl );
  // find the largest near-neighbour list, sum of square distances between members
  // as a tie-breaker
  int find_largest_nnl( std::vector<SiteCluster> &nnls );
  // put the near-neighbour lists into clusters. Destroys nnls in the process
  void cluster_nnls( std::vector<SiteCluster> &nnls ,
                     std::vector<SiteCluster> &clusters );
  // score the hpobe clusters, allowing for non-directional ones as well.
  float score_hphobe_clusters( const std::vector<SiteCluster> &clusters );
  // score the hbond clusters of type either "Donor" or "Acceptor"
  void score_hbond_clusters( std::vector<SiteCluster> &clusters ,
                             boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid ,
                             const std::string &hbond_type ,
                             std::vector<float> &scores );
  float calc_hbond_cluster_score( SiteCluster &cluster ,
                                  std::vector<SiteCluster> &virt_clusters ,
                                  boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid );
  float calc_occlusion( const SiteCluster &atom_sites ,
                        std::vector<SiteCluster> &virt_sites ,
                        boost::shared_ptr<DACLIB::VolumeGrid> &fixed_solid_grid );
  float calc_fap( float ap );
  float calc_gvp( float vp );
  // if a donor cluster has 75% of atoms in common with an acceptor cluster,
  // it's deemed a donor-acceptor cluster (probably hydroxyl groups) and only scored
  // once so remove lower-scoring one.
  void remove_don_acc_clusters( std::vector<SiteCluster> &donor_clusters ,
                                std::vector<float> &donor_scores ,
                                std::vector<SiteCluster> &acc_clusters ,
                                std::vector<float> &acc_scores ) const;

  // make new BasePPhoreSites at the positions of the virtual sites in sites
  // so we can re-use the sites clustering functions to cluster the virtual
  // sites also
  void extract_virtual_sites( const std::vector<ClusterMember> &sites ,
                              std::vector<std::vector<SinglePPhoreSite *> > &virt_sites );
  void calc_cluster_centroid( const SiteCluster &sites ,
                              double clus_centroid[3] );

  float mean_square_distance_to_centroid( const std::vector<ClusterMember> &sites );
  float mean_normal_cosine( const std::vector<ClusterMember> &sites );
  void write_scores_to_stream( std::ostream &s , const std::string &sep ) const;
  void write_scores_headers_to_stream( std::ostream &s ,
                                       const std::string &sep );
  // turn the scores into a set of string pairs, the first being the label,
  // the second being the value.
  virtual void put_scores_in_vector( std::vector<std::pair<std::string,std::string> > &scores ) const;

  std::string write_contents_to_string(); // for packing into PVM, for example.
  
  friend std::ostream &operator<<( std::ostream &s , const OverlayScore &os );

  static void set_score_method( GtplDefs::SCORE_METHOD new_method ) {
    score_method_ = new_method;
  }
  static GtplDefs::SCORE_METHOD get_score_method() { return score_method_; }

  bool operator<( const OverlayScore &rhs ) const;
  bool operator>( const OverlayScore &rhs ) const;

  // in these functions, it's not the numerical value that's less, it's what
  // that means. hence, a bigger protein clash is a poorer score, so > is used
  // in that case. In all cases, the order of molecule name is the tie-breaker,
  // so when the scores are equal the order of output is still consistent.
  bool grid_shape_tani_less( const OverlayScore &rhs ) const;
  bool gauss_shape_tanimoto_less( const OverlayScore &rhs ) const;
  bool rms_and_size_less( const OverlayScore &rhs ) const;
  bool surface_ovlp_volume_less( const OverlayScore &rhs ) const;
  bool protein_clash_less( const OverlayScore &rhs ) const;
  bool mmff_nrg_less( const OverlayScore &rhs ) const;
  // these 2 work on pareto ranking so return true if rhs dominates this
  // this one on just hbond, hphobe and volume score
  bool robins_pareto_less( const OverlayScore &rhs ) const;
  // this one on all the scores
  bool overall_pareto_less( const OverlayScore &rhs ) const;

};

// check that the two molecules and the sites are the same
bool overlay_score_names_and_sites_match( const OverlayScore *a ,
                                          const OverlayScore *b );

// so we can sort them with STL functions
class OverlayScoreMore :
    public std::binary_function<const OverlayScore * , const OverlayScore * , bool> {

public :
  
  result_type operator()( first_argument_type a , second_argument_type b ) const {
    return *a > *b;
  }

};

#endif
