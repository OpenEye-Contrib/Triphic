//
// file was_it_a_pvm_failure_message.cc
// David Cosgrove
// AstraZeneca
// 6th August 2007

#include <pvm3.h>
#include <iostream>

// Message that a PVM process has died as suggested in the PVM manual
const int TASK_DIED = 11;

namespace DACLIB {

// ********************************************************************
bool was_it_a_pvm_failure_message( int bufid , int &dead_tid ) {

  int      nbytes , msgtag;

  dead_tid = 0;
  if( -1 != pvm_bufinfo( bufid , &nbytes , &msgtag , &dead_tid ) ) {
#ifdef NOTYET
    std::cout << "bufid = " << bufid << "  nbytes = " << nbytes
              << "  msgtag = " << msgtag << "  dead_tid = " << dead_tid << std::endl;
#endif
    // the message is the tid of the dead process
    if( TASK_DIED == msgtag ) {
      pvm_upkint( &dead_tid , 1 , 1 );
      return true;
    } else {
      return false;
    }
  }

  return false;

}

}
