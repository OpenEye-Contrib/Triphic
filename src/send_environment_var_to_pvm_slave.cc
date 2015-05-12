//
// file send_environment_var_to_pvm_slave.cc
// David Cosgrove
// AstraZeneca
// 6th August 2007

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>

#include <pvm3.h>

using namespace std;

namespace DACLIB {

  // In pvm_string_subs.cc
  // pack a C++ string into a pvm buffer
  void pack_string( const string &str );

  // **************************************************************************
  void pack_environment_var( const string &var_name , const string &var_val ) {

    string msg = var_name + "=" + var_val;
    pvm_initsend( PvmDataDefault );
    pvm_pkstr( const_cast<char *>( string( "Set_Environment" ).c_str() ) );
    DACLIB::pack_string( msg );

  }

  // **************************************************************************
  // send a message to the given slave to set the given environment variable
  void send_environment_var_to_pvm_slave( const string &var_name ,
					  const string &var_val ,
					  int slave_tid ) {

    pack_environment_var( var_name , var_val );
    pvm_send( slave_tid , 0 );

  }

  // **************************************************************************
  // send a message to the given slave to set the given environment variable
  void send_environment_var_to_pvm_slaves( const string &var_name ,
                                           const string &var_val ,
                                           const vector<int> &slave_tids ) {

    pack_environment_var( var_name , var_val );
    pvm_mcast( const_cast<int *>( &slave_tids[0] ) , slave_tids.size() , 0 );

  }

  // **************************************************************************
  // receive a message about and environment variable and set it.
  void set_environment_var_from_pvm() {
    
    int      msg_len;
    char     *msg;
    
    pvm_upkint( &msg_len , 1 , 1 );
    msg = (char *) malloc( msg_len + 1 ); // putenv uses malloc/free
    pvm_upkstr( msg );
    putenv( msg );
    
  }

} // end of namespace DACLIB

