//
// file ParseSMARTSXML.H
// David Cosgrove
// AstraZeneca
// 1st May 2012
//
// This is the interface file for class ParseSMARTSXML
// which reads a SMARTS file in XML format into 2 vectors
// of pairs of strings, one holding the vector bindings,
// one holding the full definitions.

#ifndef PARSESMARTSXML_H
#define PARSESMARTSXML_H

#include <iosfwd>
#include <string>
#include <vector>

// ****************************************************************************

class ParseSMARTSXML {

public :

  ParseSMARTSXML();
  // parse_file throws a DACLIB::FileReadOpenError if it can't open the file
  void parse_file( const std::string &filename ,
                   std::vector<std::pair<std::string,std::string> > &smts ,
                   std::vector<std::pair<std::string,std::string> > &vbs );
  void parse_string( const std::string &xml ,
                     std::vector<std::pair<std::string,std::string> > &smts ,
                     std::vector<std::pair<std::string,std::string> > &vbs );
  void parse_stream( std::istream &is ,
                     std::vector<std::pair<std::string,std::string> > &smts ,
                     std::vector<std::pair<std::string,std::string> > &vbs );

  std::vector<std::pair<std::string,std::string> > smarts() const { return smarts_; }

private :

  std::vector<std::pair<std::string,std::string> > smarts_;
  std::vector<std::pair<std::string,std::string> > vbs_;

};

#endif // PARSESMARTSXML_H
