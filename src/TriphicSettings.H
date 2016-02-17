//
// file TriphicSettings.H
// David Cosgrove
// 2nd August 2007
//
// This is a public class that holds all the settings that the user used to
// start triphic with, to make them easier to pass round the place.

#ifndef DAC_TRIPHIC_SETTINGS
#define DAC_TRIPHIC_SETTINGS

#include <iosfwd>
#include <limits>
#include <string>
#include <vector>
#include <boost/program_options/options_description.hpp>

#include "gtpl_defs.H"

// **************************************************************************
class BasePPhoreSite;

class TriphicSettings {

public :

  TriphicSettings();
  // TriphicSettings( int argc , char **argv );

  // will throw a boost exception if a bad argument is given
  void parse_args( int argcc , char **argv );

  bool operator!() const;

  void print_usage( std::ostream &os ) const;

  bool set_score_method_from_primary_score() const;
  bool set_vector_usage_from_string( const std::string &str ,
                                     GtplDefs::DIRS_USAGE &du ) const;
  void set_cutoff_flags() const;

  void print_required_sites_and_points( std::ostream &os );
  // make sure the sites in req_sites_or_ and req_sites_and_ appear in query_sites
  bool check_required_sites( std::vector<BasePPhoreSite *> &query_sites );
  bool check_required_sites( std::vector<BasePPhoreSite *> &query_sites ,
                             const std::vector<std::string> &req_sites );
  // make sure the sites in req_points_or_ and req_points_and_ appear in query_sites
  bool check_required_points( std::vector<BasePPhoreSite *> &query_sites );
  bool check_required_points( std::vector<BasePPhoreSite *> &query_sites ,
                              const std::vector<std::string> &req_points );

  void send_contents_via_mpi( int dest_slave );
  void receive_contents_via_mpi(); // assumes it's from process rank 0

  std::string query_file() const { return query_file_; }
  std::vector<std::string> db_files() const { return db_files_; }
  std::string output_file() const { return output_file_; }
  std::string smarts_file() const { return smarts_file_; }
  std::string not_smarts_file() const { return not_smarts_file_; }
  std::string points_file() const { return points_file_; }
  std::string protein_file() const { return protein_file_; }
  std::string subset_file() const { return subset_file_; }
  std::string sites_score_file() const { return sites_score_file_; }
  int max_hits() const { return max_hits_; }
  bool opt_lig_rigid() const { return opt_lig_rigid_; }
  bool opt_lig_flexi() const { return opt_lig_flexi_; }
  bool gauss_shape_tanimoto() const { return gauss_shape_tani_; }
  bool grid_shape_tanimoto() const { return grid_shape_tani_; }
  bool do_omega() const { return do_omega_; }
  bool do_flipper() const { return do_flipper_; }
  bool no_warts() const { return no_warts_; }
  bool do_ions_and_tauts() const { return do_ions_and_tauts_; }
  bool all_best_hits() const { return all_best_hits_; }
  bool dont_do_sub_cliques() const { return dont_do_sub_cliques_; }
  unsigned int chunk_size() const { return chunk_size_; }
  void set_chunk_size( unsigned int ni ) { chunk_size_ = ni; }
  unsigned int max_versions() const { return max_versions_; }

  std::vector<std::string> req_sites_or() const { return req_sites_or_; }
  std::vector<std::string> req_sites_and() const { return req_sites_and_; }
  void add_req_sites_and( std::string &new_site );
  std::vector<std::string> req_points_or() const { return req_points_or_; }
  std::vector<std::string> req_points_and() const { return req_points_and_; }
  std::vector<std::string> ignore_sites() const { return ignore_sites_; }
  std::vector<std::pair<std::string,std::string> > &not_smarts_list() {
    return not_smarts_list_;
  }
  std::vector<std::string> grid_vol_files() const { return grid_vol_files_; }

  GtplDefs::SCORE_METHOD score_method() const { return score_method_; }
  GtplDefs::DIRS_USAGE  ring_norm_usage() const { return ring_norm_usage_; }
  void set_ring_norm_usage( GtplDefs::DIRS_USAGE new_du ) { ring_norm_usage_ = new_du; }
  GtplDefs::DIRS_USAGE  h_vector_usage() const { return h_vector_usage_; }
  GtplDefs::DIRS_USAGE  lp_usage() const { return lp_usage_; }
  float dist_tol() const { return dist_tol_; }
  float scaled_dist_tol() const { return scaled_dist_tol_; }
  float ring_norm_tol() const { return ring_norm_tol_; }
  float h_vector_tol() const { return h_vector_tol_; }
  float lp_tol() const { return lp_tol_; }
  float max_rms() const { return max_rms_; }
  int min_clique_size() const { return min_clique_size_; }
  float min_grid_shape_tani() const { return min_grid_shape_tani_; }
  float min_gauss_shape_tani() const { return min_gauss_shape_tani_; }
  float max_mmff_nrg() const { return max_mmff_nrg_; }
  float min_surf_vol() const { return min_surf_vol_; }
  float max_prot_clash() const { return max_prot_clash_; }
  float min_inc_vol() const { return min_inc_vol_; }
  float vol_grid_spacing() const { return vol_grid_spacing_; }
  unsigned int restart_number() const { return restart_number_; }

  bool print_sites_and_stop() const { return print_sites_and_stop_; }
  bool do_grid_shape_tani_cutoff() const { return do_grid_shape_tani_cutoff_; }
  bool do_gauss_shape_tani_cutoff() const { return do_gauss_shape_tani_cutoff_; }
  bool do_mmff_nrg_cutoff() const { return do_mmff_nrg_cutoff_; }
  bool do_surf_vol_cutoff() const { return do_surf_vol_cutoff_; }
  bool do_prot_clash_cutoff() const { return do_prot_clash_cutoff_; }
  bool do_inc_vol_cutoff() const { return do_inc_vol_cutoff_; }
  bool no_hit_conf_number() const { return no_hit_conf_number_; }
  bool output_query_to_hits() const { return output_query_to_hits_; }
  bool comma_output() const { return comma_output_; }
  bool output_scores_only() const { return output_scores_only_; }
  bool single_conf_mols() const { return single_conf_mols_; }
  bool restart() const { return restart_; }

private :

  boost::program_options::options_description desc_;

  std::string usage_text_;
  std::string query_file_;
  std::string db_files_file_; // file containing list of db files
  std::vector<std::string> db_files_;
  std::string output_file_;
  std::string smarts_file_;
  std::string not_smarts_file_; // for SMARTS that the hits mustn't match
  std::string points_file_;
  std::string protein_file_;
  std::string subset_file_;
  std::string primary_score_;
  std::string ring_norm_string_ , h_vector_string_ , lp_string_;
  std::string sites_score_file_;
  unsigned int    max_hits_;
  bool   opt_lig_rigid_ , opt_lig_flexi_;
  bool   gauss_shape_tani_;
  mutable bool grid_shape_tani_;
  bool   do_omega_ , do_flipper_;
  bool   do_ions_and_tauts_;
  bool   no_warts_;
  bool   all_best_hits_; /* keep the best-scoring hit for each clique, rather
                            than the best-scoring hit for each target. Doesn't
                            over-ride max_hits_. */
  bool dont_do_sub_cliques_; /* when doing a clique, don't also generate all
                                sub-cliques of at least size min_clique_size_.
                                Default is to generate them. */
  bool restart_; // whether it's a restart job
  unsigned int chunk_size_; // for parallel jobs
  unsigned int max_versions_; // for a particular molecule, max number of ions, tautomers, stereoisomers

  std::vector<std::string> req_sites_or_ , req_sites_and_;
  std::vector<std::string> req_points_or_ , req_points_and_;
  std::vector<std::string> ignore_sites_;
  std::vector<std::pair<std::string,std::string> > not_smarts_list_;
  std::vector<std::string> grid_vol_files_;

  mutable GtplDefs::SCORE_METHOD score_method_;
  mutable GtplDefs::DIRS_USAGE  ring_norm_usage_ , h_vector_usage_ , lp_usage_;
  float dist_tol_; // tolerance for matching dists in cliques
  float scaled_dist_tol_; // tolerance will be scaled_dist_tol_ * query distance.
  float ring_norm_tol_ , h_vector_tol_ , lp_tol_; /* tolerances for vectors */
  float max_rms_; // max rms allowed after overlay
  int   min_clique_size_;
  float min_grid_shape_tani_ , min_gauss_shape_tani_;
  float max_mmff_nrg_ , min_surf_vol_ , max_prot_clash_;
  float min_inc_vol_;
  float vol_grid_spacing_;
  unsigned int restart_number_; // molecule number for restart job
  mutable bool do_grid_shape_tani_cutoff_ , do_gauss_shape_tani_cutoff_;
  mutable bool do_mmff_nrg_cutoff_ , do_surf_vol_cutoff_;
  mutable bool do_prot_clash_cutoff_ , do_inc_vol_cutoff_;
  bool no_hit_conf_number_ , output_query_to_hits_ , print_sites_and_stop_;
  bool comma_output_ , output_scores_only_ , single_conf_mols_;

  void build_program_options();
  void read_database_files_file();

};

#endif
