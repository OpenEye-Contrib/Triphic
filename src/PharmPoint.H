//
// file PharmPoint.H
// David Cosgrove
// AstraZeneca
// 12th November 2002
//
// This is the interface file for the class PharmPoint.  It holds definitions
// for pharmacophore points used by various programsincluding Triphic, Loob
// and Plurality.

#ifndef DAC_PHARM_POINT__
#define DAC_PHARM_POINT__

#include <map>
#include <set>
#include <string>
#include <vector>

namespace OEChem {
  class OESubSearch;
}

// ******************************************************************************

class PharmPoint {

public:

  PharmPoint() { set_defaults(); }

  virtual ~PharmPoint() {}

  std::map<std::string,std::vector<std::string> > &points_defs() { return points_defs_; }
  std::vector<std::string> &h_vector_points() { return h_vector_points_; }
  std::vector<std::string> &lp_vector_points() { return lp_vector_points_; }
  std::vector<std::string> &ring_normal_points() { return ring_normal_points_; }
  std::set<std::string> &unique_smarts() { return unique_smarts_; }
  std::vector<int> &smarts_to_points() { return smarts_to_points_; }
  bool &hphobes_itmoc() { return hphobes_itmoc_; }
  bool &hphobes_itmoc_alo() { return hphobes_itmoc_alo_; }
  std::string &itmoc_label() { return itmoc_label_; }
  std::string &itmoc_alo_label() { return itmoc_alo_label_; }

  // file readers throw DACLIB::FileReadOpenError if they feel the urge
  virtual void read_points_file( const std::string &file_name );
  virtual void read_points_xml_file( const std::string &file_name );
  virtual void read_points_xml_string( const std::string xml );
  virtual void read_points_xml_stream( std::istream &is );

  virtual void report_points_defined( std::ostream &os ) const;

  virtual void clear_data();

  // make sure all the SMARTS defs in the points_defs are in the input_smarts.
  // throws a DACLIB::SMARTSDefnError listing the first missing one if there
  // is one
  virtual void check_points_smarts( const std::vector<std::pair<std::string,std::string> > &input_smarts );

  // take the string, a points type name, and return the corresponding integer
  // code
  int type_code_from_string( const std::string &type_string ) const;

protected:

  std::map<std::string,std::vector<std::string> > points_defs_;
  std::vector<std::string> h_vector_points_;
  std::vector<std::string> lp_vector_points_;
  std::vector<std::string> ring_normal_points_;
  std::set<std::string> unique_smarts_;
  std::vector<int> smarts_to_points_; // mapping of the points that each smarts is in
  bool        hphobes_itmoc_; /* in the manner of Catalyst (JCICS,34,1297-1308,
				 (1994)) */
  bool        hphobes_itmoc_alo_; // aliphatics only
  std::string      itmoc_label_ , itmoc_alo_label_;

  virtual void set_defaults();

};

void build_oesubsearches( PharmPoint &pharm_points ,
			  const std::vector<std::pair<std::string,std::string> > &smarts ,
			  std::map<std::string,OEChem::OESubSearch *> &oesubs );

#endif
