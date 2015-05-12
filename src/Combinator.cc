//
// file Combinator.cc
// David Cosgrove
// AstraZeneca
// 11th July 2008

#include "Combinator.H"

#include <iostream>
#include <iterator>

using namespace std;

// ********************************************************************
void Combinator::reset() {

  for( int i = 0 ; i < k_ ; ++i ) {
    curr_comb_.push_back( i );
  }

}

// ********************************************************************
void Combinator::step() {

  if( at_end() ) {
    return;
  }

  // if k_ is zero, we shouldn't get here so k_ - 1 is safe.
  for( int i = k_ - 1 ; i >= 0 ; --i ) {
    if( curr_comb_[i] < n_ - k_ + i ) {
      for( int j = curr_comb_[i] ; i < k_ ; ++i ) {
	curr_comb_[i] = ++j;
      }
      return;
    }
  }

  // if we're here, we can't step further
  curr_comb_.clear();

}
