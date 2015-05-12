//
// file PluralitySettings.H
// Dave Cosgrove
// 13th September 2007
//
// This is a public class that holds all the settings that the user used to
// start plurality with, to make them easier to pass round the place.

#ifndef DAC_PLURALITY_SETTINGS
#define DAC_PLURALITY_SETTINGS

#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

#include <boost/program_options/options_description.hpp>

#include "gtpl_defs.H"

typedef enum { BEST_HITS_ONLY , ONE_HIT_ONLY , ALL_HITS } HITS_TO_OUTPUT;

// **************************************************************************
class PluralitySettings {

public :

  PluralitySettings();
  PluralitySettings( int argc , char **argv );

  bool operator!() const;
  void print_usage( std::ostream &os ) const;

  void pack_contents_into_pvm_buffer();
  void unpack_contents_from_pvm_buffer();

  std::string pphore_file() const { return pphore_file_; }
  std::vector<std::string> db_files() const { return db_files_; }
  std::string output_file() const { return output_file_; }
  std::string smarts_file() const { return smarts_file_; }
  std::string points_file() const { return points_file_; }
  std::string protein_file() const { return protein_file_; }
  std::string subset_file() const { return subset_file_; }
  std::string not_smarts_file() const { return not_smarts_file_; }

  std::vector<std::pair<std::string,std::string> > &not_smarts_list() { return not_smarts_list_; }
  std::vector<std::string> &grid_vol_files() { return grid_vol_files_; }
  bool scores_only() const { return scores_only_; }
  bool comma_output() const { return comma_output_; }
  HITS_TO_OUTPUT hits_to_output() const { return hits_to_output_; }

  std::string slave_name() const { return slave_name_; }
  int num_slave_procs() const { return num_slave_procs_; }
  std:: string pvm_hosts_file() const { return pvm_hosts_file_; }

private :

  std::string pphore_file_;
  std::string db_files_file_; // file containing list of db files
  std::vector<std::string> db_files_;
  std::string output_file_;
  std::string smarts_file_;
  std::string points_file_;
  std::string protein_file_;
  std::string subset_file_;
  std::string not_smarts_file_;
  std::vector<std::string> grid_vol_files_;

  std::vector<std::pair<std::string,std::string> > not_smarts_list_;
  bool scores_only_;
  bool comma_output_;
  std::string hits_to_output_str_;
  mutable HITS_TO_OUTPUT hits_to_output_;

  // for parallel runs
  std::string slave_name_;
  int    num_slave_procs_;
  std::string pvm_hosts_file_;

  std::string usage_text_;

  bool set_hits_to_output_from_string() const;
  void build_program_options( boost::program_options::options_description &desc );
  void read_database_files_file();

};

#endif
