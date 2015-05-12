//
// file LoobSettings.H
// Dave Cosgrove
// AstraZeneca
// 2nd October 2007
//
// This is a public class that holds all the settings used to start loob.

#ifndef DAC_LOOB_SETTINGS
#define DAC_LOOB_SETTINGS

#include <string>
#include <vector>
#include <boost/program_options/options_description.hpp>

// *******************************************************************

class LoobSettings {

public :

  LoobSettings( int argc , char **argv );
  
  std::string usage_text_;
  std::string db_file_;
  std::string ascii_fp_file_ , compact_fp_file_ , names_fp_file_;
  std::string log_file_ , bit_names_file_;
  std::string smarts_file_;
  std::string points_file_;
  std::string subset_file_;
  std::string bit_separator_;
  bool   pairs_ , triplets_ , quadruplets_;
  bool   single_conf_mols_;
  bool   dont_compress_fps_;
  bool   chiral_fps_;
  bool   ascii_fps_;
  bool   compact_fps_;
  bool   labels_not_bits_;
  bool   bit_labels_to_ascii_file_;
  int    compressed_counts_cutoff_;

  std::vector<float> dist_bounds_ , sq_dist_bounds_;

private :

  void build_program_options( boost::program_options::options_description &desc );
  void polish_args();

};

#endif
