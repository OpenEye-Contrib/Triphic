//
// file find_cliques.cc
// David Cosgrove
// AstraZeneca
// 10th April 2006
//
// This file contains the functions for finding cliques in two vectors of
// PPhoreSites, containing sites of the types given.
// If scaled_dist_tol is -1.0F, uses the dist_tol value to decide if a distance
// between two query sites and two target sites is a match.  If not, then it
// takes the tolerance to be scaled_dist_tol_ * query_dist.

#include <algorithm>
#include <iostream>

#include "stddefs.H"
#include "Combinator.H"
#include "SinglePPhoreSite.H"

using namespace std;

extern "C" {

void bron_kerbosch( char **adj_table , int num_nodes ,
                    unsigned short ***cliques , int **clique_sizes ,
                    int *num_cliques , int *min_clique_size );

};

class SortClique : public binary_function<const pair<int,int> &,
    const pair<int,int> &, bool> {
public :
  result_type operator()( first_argument_type a, second_argument_type b ) const {
    return a.first < b.first;
  }

};

class SortCliques : public binary_function<const vector<int> &,
    const vector<int> &, bool> {
public :
  result_type operator()( first_argument_type a, second_argument_type b ) const {

    if( a.size() == b.size() ) {
      for( int i = 0 , is = a.size() ; i < is ; i += 2 ) {
        if( a[i] != b[i] ) {
          return a[i] < b[i];
        }
      }
      // we're here, the first set are all the same, go on the second
      for( int i = 1 , is = a.size() ; i < is ; i += 2 ) {
        if( a[i] != b[i] ) {
          return a[i] < b[i];
        }
      }
      // we're here, they must be the same which means a !< b
      return false;
    } else {
      return a.size() < b.size();
    }

  }

};

// ****************************************************************************
void find_matching_pairs( const vector<BasePPhoreSite *> &sites1 ,
                          const vector<SinglePPhoreSite *> &sites2 ,
                          const vector<int> &allowed_types , int &num_pairs ,
                          int *&pairs ) {

  pairs = new int[2 * sites1.size() * sites2.size()];
  num_pairs = 0;
  for( unsigned int i = 0 ; i < allowed_types.size() ; ++i ) {
    for( unsigned int j = 0 ; j < sites1.size() ; ++j ) {
      if( sites1[j]->get_type_code() != allowed_types[i] ) {
        continue;
      }
      for( unsigned int k = 0 ; k < sites2.size() ; ++k ) {
        if( sites2[k]->get_type_code() == allowed_types[i] ) {
          pairs[2 * num_pairs] = j;
          pairs[2 * num_pairs + 1] = k;
          ++num_pairs;
        }
      }
    }
  }

}

// ****************************************************************************
void find_matching_pairs( const vector<BasePPhoreSite *> &sites1 ,
                          const vector<SinglePPhoreSite *> &sites2 ,
                          int &num_pairs , int *&pairs ) {

  pairs = new int[2 * sites1.size() * sites2.size()];
  num_pairs = 0;
  for( unsigned int j = 0 ; j < sites1.size() ; ++j ) {
    for( unsigned int k = 0 ; k < sites2.size() ; ++k ) {
      if( sites1[j]->get_type_code() == sites2[k]->get_type_code() ) {
        pairs[2 * num_pairs] = j;
        pairs[2 * num_pairs + 1] = k;
        ++num_pairs;
      }
    }
  }

}

// ****************************************************************************
void build_adj_table( const vector<BasePPhoreSite *> &sites1 ,
                      const vector<SinglePPhoreSite *> &sites2 ,
                      int &num_pairs , int *&pairs ,
                      float dist_tol , float scaled_dist_tol ,
                      char **&adj_table ) {

  bool scaled_tol = scaled_dist_tol < 0.0F ? false : true;

  DACLIB::make_square_matrix( adj_table , num_pairs );
  fill( adj_table[0] , adj_table[0] + num_pairs * num_pairs , 0 );
  for( int i = 0 ; i < num_pairs - 1 ; ++i ) {
    adj_table[i][i] = 1; // BK algorithm needs diagonal set to 1
    int *pair_i = pairs + 2 * i;
    for( int j = i + 1 ; j < num_pairs ; ++j ) {
      int *pair_j = pairs + 2 * j;
      if( pair_i[0] == pair_j[0] || pair_i[1] == pair_j[1] ) {
        continue; // can't have point matching itself
      }
      float dist1 = sites1[pair_i[0]]->distance( sites1[pair_j[0]]->coords() );
      float dist2 = sites2[pair_i[1]]->distance( sites2[pair_j[1]]->coords() );
      if( scaled_tol ) {
        dist_tol = scaled_dist_tol * dist1;
      }
      if( fabs( dist1 - dist2 ) < dist_tol ) {
        adj_table[i][j] = adj_table[j][i] = 1;
      }
    }
  }
  adj_table[num_pairs-1][num_pairs-1] = 1;

}

// ****************************************************************************
// sort the clique into ascending order of first pair number
void sort_clique( vector<int> &clique ) {

  vector<pair<int,int> > idx;
  for( unsigned int i = 0 , is = clique.size() ; i < is ; i += 2 ) {
    idx.push_back( make_pair( clique[i] , i/2 ) );
  }
  sort( idx.begin() , idx.end() , SortClique() );

  vector<int> out_clique;
  for( unsigned int i = 0 , is = idx.size() ; i < is ; ++i ) {
    out_clique.push_back( clique[2*idx[i].second] );
    out_clique.push_back( clique[2*idx[i].second+1] );
  }

  clique = out_clique;

}

// ****************************************************************************
// for cliques that are greater in size than min_clique_size, we need to
// create all the sub-cliques. These won't score better than the clique on
// Size/RMS, but they might pass a filter that the whole clique doesn't, or
// we might be using a different primary score where the sub-clique gives a
// better scoring overlay.
void create_sub_cliques( int min_clique_size ,
                         vector<vector<int> > &cliques ) {

  unsigned int twice_mcs = 2 * min_clique_size;
  vector<vector<int> > new_cliques;
  for( int i = 0 , is = cliques.size() ; i < is ; ++i ) {
    if( cliques[i].size() > twice_mcs ) {
      for( int k = min_clique_size , ks = cliques[i].size() / 2 ; k < ks ; ++k ) {
        Combinator comb( k , ks );
        while( !comb.at_end() ) {
          const vector<int> &next_comb = comb.value();
          cliques.push_back( vector<int>() );
          for( int j = 0 , js = next_comb.size() ; j < js ; ++j ) {
            cliques.back().push_back( cliques[i][2 * next_comb[j]] );
            cliques.back().push_back( cliques[i][2 * next_comb[j] + 1] );
          }
          comb.step();
        }
      }
    }
  }

  sort( cliques.begin() , cliques.end() , SortCliques() );
  // different cliques might give rise to the same sub-clique. It's clearly
  // wasted effort scoring them more than once each.
  cliques.erase( unique( cliques.begin() , cliques.end() ) , cliques.end() );

}

// ****************************************************************************
void find_cliques( const vector<BasePPhoreSite *> &sites1 ,
                   const vector<SinglePPhoreSite *> &sites2 ,
                   float dist_tol , float scaled_dist_tol ,
                   bool dont_do_sub_cliques ,
                   int min_clique_size , vector<vector<int> > &cliques ) {

  // find all matching pairs
  int num_pairs , *pairs = 0;
  find_matching_pairs( sites1 , sites2 , num_pairs , pairs );

  if( !num_pairs ) {
    delete [] pairs;
    return;
  }

  // build the adjacency table for the pairs
  char **adj_table = 0;
  build_adj_table( sites1 , sites2 , num_pairs , pairs , dist_tol ,
                   scaled_dist_tol , adj_table );

  // find the cliques
  int      num_cliques , *clique_sizes;
  unsigned short **raw_cliques;
  bron_kerbosch( adj_table , num_pairs , &raw_cliques , &clique_sizes ,
                 &num_cliques , &min_clique_size );

  // translate the raw cliques to appropriate return types
  cliques.clear();
  for( int i = 0 ; i < num_cliques ; ++i ) {
    if( clique_sizes[i] >= min_clique_size ) {
      cliques.push_back( vector<int>() );
      for( int j = 0 ; j < clique_sizes[i] ; ++j ) {
        cliques.back().push_back( pairs[2 * raw_cliques[i][j]] );
        cliques.back().push_back( pairs[2 * raw_cliques[i][j] + 1] );
      }
      sort_clique( cliques.back() );
    }
    free( raw_cliques[i] ); // b-k is C code, uses malloc/calloc
  }
  free( raw_cliques );
  free( clique_sizes );

  sort( cliques.begin() , cliques.end() , SortCliques() );

  if( !dont_do_sub_cliques ) {
    create_sub_cliques( min_clique_size , cliques );
  }

  delete [] pairs;
  delete [] adj_table[0];
  delete [] adj_table;

}
