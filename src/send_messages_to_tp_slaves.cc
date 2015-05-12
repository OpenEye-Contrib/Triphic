//
// file send_messages_to_tp_slaves.cc
// David Cosgrove
// 25th September 2007
// Various functions for sending messages to the triphic and plurality slaves

#include <iostream>
#include <string>
#include <vector>

#include <pvm3.h>

#include <sys/param.h>
#include <unistd.h>

#include <oechem.h>

using namespace std;
using namespace OEChem;

namespace DACLIB {
bool was_it_a_pvm_failure_message( int bufid , int &dead_tid );
void send_environment_var_to_pvm_slaves( const string &var_name ,
                                         const string &var_val ,
                                         const vector<int> &slave_tids );
void pack_string( const string &str );
void unpack_string( string &str );
}

// in eponymous file
bool step_oemolstream( oemolistream &ims , int step );

// ********************************************************************
// spell from t'Interweb
string getcwd(){

  char *buffer = new char[MAXPATHLEN];
  getcwd(buffer,MAXPATHLEN);
  if(buffer != NULL){
    string ret(buffer);
    delete[] buffer;
    return ret;
  } else {
    return string();
  }

}

// ********************************************************************
void get_done_message_from_slave( vector<int> &slave_tids ,
                                  vector<char> &slave_busy ) {

  int bufid = pvm_recv( -1 , -1 );
  int done_tid;
  if( DACLIB::was_it_a_pvm_failure_message( bufid , done_tid ) ) {

    cerr << "Process " << done_tid << " has gone belly up, taking all" << endl
         << "its results with it. Carrying on, but there will be missing"
         << endl
         << "hits." << endl;
    cout << "Process " << done_tid << " has gone belly up, taking all" << endl
         << "its results with it. Carrying on, but there will be missing"
         << endl
         << "hits." << endl;
    for( int i = 0 , is = slave_tids.size() ; i < is ; ++i ) {
      if( done_tid == slave_tids[i] ) {
        slave_tids.erase( slave_tids.begin() + i );
        slave_busy.erase( slave_busy.begin() + i );
        break;
      }
    }
    if( slave_tids.empty() ) {
      cerr << "All slaves are now dead, so that's it for now. This could" << endl
           << "be a problem with the program, it could be a problem with" << endl
           << "your database, or it could be a problem with the machines" << endl
           << "the slaves were running on, including the possibility that" << endl
           << "the OpenEye license isn't set correctly on the slave" << endl
           << "machines." << endl;
      pvm_exit();
      exit( 1 );
    }
    // still need to get another slave that's done.
    get_done_message_from_slave( slave_tids , slave_busy );

  } else {

    // if we're here, the message wasn't a failure message. It could be
    // a progress report of a done message. done_tid isn't
    // filled in this case, so need to get it from the message, flag the slave
    // as free and we're done.
    char msg[1000]; // it'll be big enough for the message header
    pvm_upkstr( msg );
    cout << "Message : " << msg << endl;
    if( !strcmp( msg , "Finished Searching" ) ) {
      pvm_upkint( &done_tid , 1 , 1 );
      for( int i = 0 , is = slave_tids.size() ; i < is ; ++i ) {
        if( done_tid == slave_tids[i] ) {
          slave_busy[i] = 0;
          break;
        }
      }
      cout << "Number of slaves still to report : "
           << count( slave_busy.begin() , slave_busy.end() , 1 ) << endl;
    } else if( !strcmp( "Progress Report" , msg ) ) {
      pvm_upkint( &done_tid , 1 , 1 );
      string prog_rep;
      DACLIB::unpack_string( prog_rep );
      cout << done_tid << " : " << prog_rep;
    }

  }

}

// ********************************************************************
void pack_smarts_defs_into_pvm_buffer( vector<pair<string,string> > &input_smarts ,
                                       vector<pair<string,string> > &smarts_sub_defn ) {

  int num_to_send = input_smarts.size();
  pvm_pkint( &num_to_send , 1 , 1 );
  for( int i = 0 ; i < num_to_send ; ++i ) {
    DACLIB::pack_string( input_smarts[i].first );
    DACLIB::pack_string( input_smarts[i].second );
  }

  num_to_send = smarts_sub_defn.size();
  pvm_pkint( &num_to_send , 1 , 1 );
  for( int i = 0 ; i < num_to_send ; ++i ) {
    DACLIB::pack_string( smarts_sub_defn[i].first );
    DACLIB::pack_string( smarts_sub_defn[i].second );
  }

}

// ********************************************************************
void unpack_smarts_defs_from_pvm_buffer( vector<pair<string,string> > &input_smarts ,
                                         vector<pair<string,string> > &smarts_sub_defn ) {

  int num_to_rec;

  pvm_upkint( &num_to_rec , 1 , 1 );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    input_smarts.push_back( make_pair( string( "" ) , string( "" ) ) );
    DACLIB::unpack_string( input_smarts[i].first );
    DACLIB::unpack_string( input_smarts[i].second );
  }

  pvm_upkint( &num_to_rec , 1 , 1 );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    smarts_sub_defn.push_back( make_pair( string( "" ) , string( "" ) ) );
    DACLIB::unpack_string( smarts_sub_defn[i].first );
    DACLIB::unpack_string( smarts_sub_defn[i].second );
  }

}

// ********************************************************************
void send_openeye_license_to_slaves( const vector<int> &slave_tids ) {

  if( getenv( "OE_LICENSE" ) ) {
    string oe_lic( getenv( "OE_LICENSE" ) );
    DACLIB::send_environment_var_to_pvm_slaves( string( "OE_LICENSE" ) ,
                                                oe_lic , slave_tids );
  }

}

// ********************************************************************
void send_cwd_to_slaves( const vector<int> &slave_tids ) {

  string cwd = getcwd();
  if( cwd.length() ) {
    pvm_initsend( PvmDataDefault );
    pvm_pkstr( const_cast<char *>( string( "New_CWD" ).c_str() ) );
    DACLIB::pack_string( cwd );
    pvm_mcast( const_cast<int *>( &slave_tids[0] ) , slave_tids.size() , 0 );
  }

}

// ********************************************************************
void send_next_targets_to_slave( const vector<vector<unsigned char> > &mol_recs ,
                                 int slave_tid ) {

  pvm_initsend( PvmDataDefault );
  pvm_pkstr( const_cast<char *>( string( "Next_Targets" ).c_str() ) );
  int num_to_send = mol_recs.size();
  pvm_pkint( &num_to_send , 1 , 1 );
  for( int i = 0 ; i < num_to_send ; ++i ) {
    int num_bytes_to_send = mol_recs[i].size();
    pvm_pkint( &num_bytes_to_send , 1 , 1 );
    pvm_pkbyte( (char *)( &mol_recs[i][0] ) , num_bytes_to_send , 1 );
  }

  pvm_send( slave_tid , 0 );

}

// ********************************************************************
void send_next_targets_to_slave( const vector<vector<unsigned char> > &mol_recs ,
                                 vector<int> &slave_tids ,
                                 vector<char> &slave_busy ) {

  // find a slave to use, if we do, send the records and return
  for( int i = 0 , is = slave_busy.size() ; i < is ; ++i ) {
    if( !slave_busy[i] ) {
      send_next_targets_to_slave( mol_recs , slave_tids[i] );
      slave_busy[i] = 1;
      return;
    }
  }

  // if we're here, all slaves are busy. Wait for one to come free, then
  // re-send. get_done_message_from_slave may change slave_tids and slave_busy
  // if the slave has died. If all the slaves die, the program aborts.
  get_done_message_from_slave( slave_tids , slave_busy );
  send_next_targets_to_slave( mol_recs , slave_tids , slave_busy );

  if( slave_tids.size() < 1 )
    exit( 1 );

}

// ********************************************************************
void send_finished_messages( const vector<int> &slave_tids ) {

  pvm_initsend( PvmDataDefault );
  pvm_pkstr( const_cast<char *>( string( "Finished" ).c_str() ) );
  pvm_mcast( const_cast<int *>( &slave_tids[0] ) , slave_tids.size() , 0 );

}

// ********************************************************************
void send_done_message_to_master( int master_tid , int my_tid ) {

  pvm_initsend( PvmDataDefault );
  pvm_pkint( &my_tid , 1 , 1 );
  pvm_pkstr( const_cast<char *>( string( "Finished Searching" ).c_str() ) );
  pvm_send( master_tid , 0 );

}

// ********************************************************************
void parse_oebin_records( const vector<vector<unsigned char> > &mol_recs ,
                          vector<OEMol *> &mols ) {

  for( int i = 0 , is = mol_recs.size() ; i < is ; ++i ) {
    if( mol_recs[i].empty() )
      break;
    oemolistream ims;
    ims.SetFormat( OEFormat::OEB );
    ims.openstring( &mol_recs[i][0] , mol_recs[i].size() );
    mols.push_back( new OEMol );
    ims >> *mols[i];
    OEAssignAromaticFlags( *mols[i] , OEAroModelDaylight );
  }

}

// ********************************************************************
void receive_next_target( vector<OEMol *> &target_mols ) {

  int num_to_rec;
  pvm_upkint( &num_to_rec , 1 , 1 );
  vector<vector<unsigned char> > mol_recs;
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    int num_bytes_to_rec;
    pvm_upkint( &num_bytes_to_rec , 1 , 1 );
    mol_recs.push_back( vector<unsigned char>( num_bytes_to_rec , 0 ) );
    pvm_upkbyte( (char *) &mol_recs[i][0] , num_bytes_to_rec , 1 );
  }

  parse_oebin_records( mol_recs , target_mols );

}

// ********************************************************************
void send_database_details_to_slaves( const vector<int> &slave_tids ) {

  static const string db_steps_msg( "Database_Steps" );
  for( int i = 0 , is = slave_tids.size() ; i < is ; ++i ) {
    pvm_initsend( PvmDataDefault );
    pvm_pkstr( const_cast<char *>( db_steps_msg.c_str() ) );
    // send the start molecule for this slave, and the step size
    pvm_pkint( &i , 1 , 1 );
    pvm_pkint( &is , 1 , 1 );
    pvm_send( slave_tids[i] , 0 );
  }

}

// ********************************************************************
void receive_new_cwd() {

  string new_cwd;
  DACLIB::unpack_string( new_cwd );
  chdir( new_cwd.c_str() );

}

// ********************************************************************
void receive_database_details( int &db_start , int &db_step ) {

  pvm_upkint( &db_start , 1 , 1 );
  pvm_upkint( &db_step , 1 , 1 );

}

// **************************************************************************
void send_progress_to_master( const string &progress_report ) {

  pvm_initsend( PvmDataDefault );
  const static string msg_header( "Progress Report" );
  pvm_pkstr( const_cast<char *>( msg_header.c_str() ) );
  int my_tid = pvm_mytid() , master_tid = pvm_parent();
  pvm_pkint( &my_tid , 1 , 1 );
  DACLIB::pack_string( progress_report );
  pvm_send( master_tid , 0 );

}
