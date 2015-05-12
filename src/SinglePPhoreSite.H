//
// file SinglePPhoreSite.H
// David Cosgrove
// AstraZeneca
// 10th October 2006
//
// This is the interface for class SinglePPhoreSite, derived from BasePPhoreSite.
// It uses 1 or more atoms from a single molecule to make a pharmacophore site.

#ifndef DAC_SINGLE_PPHORE_SITE__
#define DAC_SINGLE_PPHORE_SITE__

#include "BasePPhoreSite.H"

#include <boost/shared_ptr.hpp>

// maximum number of atoms that can be in a site. If it's larger than this,
// we'll have to split it up. Needs to be able to cope with a benzene ring with
// all heavy-atom substituents.
static const int MAX_SITE_ATOMS = 12;

// ****************************************************************************

class SinglePPhoreSite : public BasePPhoreSite {

protected :

  std::string   parent_mol_name_; // name of molecule site is derived from

  // the atoms that were used to define the site - contains the OEChem index
  // i.e. the results of GetIdx() on the constituent atoms if the sites were
  // built from OEChem molecules. This is more often than not the same as the
  // sequential order number of the atom in the molecule, but it's not
  // guaranteed.
  int      num_site_atoms_;
  int      site_atoms_[MAX_SITE_ATOMS];

  void reset_data();
  void copy_data( const SinglePPhoreSite &c );

public :

  SinglePPhoreSite();
  // tc and ts are the type code and type string respectively. Throws an int
  // exception if num_atoms is too large, giving the maximum number
  SinglePPhoreSite( const double cds[3] , const double dir[3] , int tc , const std::string &ts ,
  const std::string &lbl , bool use_dir , int num_atoms ,
  const int *atoms , const std::string &parent_mol_name );
  // this version builds the SinglePPhoreSite without the atoms and therefore won't
  // throw an exception.
  SinglePPhoreSite( const double cds[3] , const double dir[3] , int tc ,
  const std::string &ts ,
  const std::string &lbl , bool use_dir ,
  const std::string &parent_mol_name );
  SinglePPhoreSite( const SinglePPhoreSite &c ) { copy_data( c ); };

  virtual ~SinglePPhoreSite() {}

  SinglePPhoreSite &operator=( const SinglePPhoreSite &c ) {
    if( this == &c )
      return *this;
    copy_data( c );
    return *this;
  }
  
  virtual std::string get_full_name() const;
  std::string parent_mol_name() const { return parent_mol_name_; }
  int num_site_atoms() const { return num_site_atoms_; }
  const int *site_atoms() const { return site_atoms_; }

  void write_to_stream( std::ostream &os ) const;
  void read_from_stream( std::istream &is );
  virtual void brief_report( std::ostream &os ) const;

  virtual std::string type_string_for_stream() {
    return std::string( "Single_PPhore_Site" );
  }

  friend std::ostream &operator<<( std::ostream &os ,
                                   const SinglePPhoreSite &site );
  friend std::istream &operator>>( std::ostream &is ,
                                   SinglePPhoreSite &site );

};

float calc_overlay_trans( const std::vector<BasePPhoreSite *> &sites1 , 
                          const std::vector<SinglePPhoreSite *> &sites2 ,
                          const std::vector<int> &pairs ,
                          OverlayTrans &overlay_trans , bool use_ring_norm_dirs ,
                          bool use_h_vec_dirs , bool use_lp_dirs );
float calc_overlay_trans( const std::vector<SinglePPhoreSite *> &sites1 , 
                          const std::vector<SinglePPhoreSite *> &sites2 ,
                          const std::vector<int> &pairs ,
                          OverlayTrans &overlay_trans , bool use_ring_norm_dirs ,
                          bool use_h_vec_dirs , bool use_lp_dirs );

float calculate_site_rms( const std::vector<BasePPhoreSite *> &sites1 ,
                          const std::vector<SinglePPhoreSite *> &sites2 ,
                          const std::vector<int> &pairs );
float calculate_site_rms( const std::vector<BasePPhoreSite *> &sites1 ,
                          const std::vector<boost::shared_ptr<SinglePPhoreSite> > &sites2 ,
                          const std::vector<int> &pairs );

void overlay( SinglePPhoreSite &site , const OverlayTrans &ot );
void overlay( std::vector<SinglePPhoreSite *> &sites , const OverlayTrans &ot );
void overlay( std::vector<boost::shared_ptr<SinglePPhoreSite> > &sites , const OverlayTrans &ot );

// confirm that the sites named in the given clique have directions within the
// tolorance. The sites will already have been overlaid. The target sites may
// be twiddled if this is possible, to attempt to line them up with the query.
bool confirm_site_vectors( const std::vector<int> &clique ,
                           const std::vector<BasePPhoreSite *> &query_sites ,
                           std::vector<boost::shared_ptr<SinglePPhoreSite> > &target_sites ,
                           bool check_ring_norms , float ring_norm_tol ,
                           bool check_h_vectors , float h_vector_tol ,
                           bool check_lps , float lp_tol );
bool confirm_ring_norms( const BasePPhoreSite &site1 ,
                         SinglePPhoreSite &site2 , float tol );
bool confirm_dir_vectors( const BasePPhoreSite &site1 ,
                          SinglePPhoreSite &site2 , float tol ,
                          GtplDefs::DIRS_TYPE dirs_type );

#endif
