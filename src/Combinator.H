//
// file Combinator.H
// David Cosgrove
// AstraZeneca
// 11th July 2008
//
// This class steps through all combinations of k things
// selected from n. It is ripped out of Flamig's
// 'Practical Algorithms in C++' pp102 but I would like to think that this
// is a product of laziness rather than inability!
// It doesn't check whether k and n are passed in as valid, for example if
// k is greater than n. Perhaps it should throw an exception, but for now
// it gives k from k in that case.

#ifndef DAC_COMBINATOR
#define DAC_COMBINATOR

#include <vector>

// **********************************************************************

class Combinator {

public :

  // default - choose nothing from nothing
  Combinator() : k_( 0 ) , n_( 0 ) {}
  // select k things from n. k must be <= n, but we won't be checking.
  Combinator( int k , int n ) : k_( k ) , n_( n ) {
    reset();
  }
  ~Combinator() {}

  void reset();
  const std::vector<int> &value() const { return curr_comb_; }
  void step();
  bool at_end() const {
    // when we step beyond the end, we'll empty curr_comb_
    return curr_comb_.empty();
  }
  int num_to_choose() const { return k_; }
  int num_to_choose_from() const { return n_; }

private :

  int k_; // number to choose
  int n_; // number to choose from
  std::vector<int> curr_comb_;

};

#endif
