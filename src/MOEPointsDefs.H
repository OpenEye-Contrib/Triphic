//
// file MOEPointsDefs.H
// David Cosgrove
// AstraZeneca
// 23rd March 2015
//
// This file contains a set of points definitions for
// triphic for use with MOE pharmacophore files.
// They're clearly not exact matches to the MOE definitions,
// but hopefully very similar. They are re-labellings of the ones
// in DefaultPointsDefs.H
//

#ifndef MOEPOINTSDEFS_H
#define MOEPOINTSDEFS_H

#include <string>

namespace DACLIB {

const std::string MOE_POINT_DEFS = "\n\
<?xml version=\"1.0\"?>\n\
<features>\n\
  <!-- vector bindings used in pharmacophore features -->\n\
  <smarts>\n\
    <name>Ylid</name>\n\
    <value>[A-]~[N,n;+;H0]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>AllDon</name>\n\
    <value>[O,N,n;!H0]</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ac1</name>\n\
    <value>[O;X1]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ac2</name>\n\
    <value>[OH]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ac3</name>\n\
    <value>O([CX4])[CX4]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>OAc</name>\n\
    <value>[$Ac1,$Ac2,$Ac3]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ac4</name>\n\
    <value>[nX2]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ac5</name>\n\
    <value>[NX2]=C</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ac6</name>\n\
    <value>[N;X3;H2][CX4]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ac7</name>\n\
    <value>[N;X3;H]([CX4])[CX4]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ac8</name>\n\
    <value>[N;X3]([CX4])([CX4])[CX4]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ac9</name>\n\
    <value>[N;X1]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>NonAccN</name>\n\
    <value>[$(ncC=O)]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>NAc</name>\n\
    <value>[$Ac4,$Ac5,$Ac6,$Ac7,$Ac8,$Ac9;!$NonAccN]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ne1</name>\n\
    <value>[*;-]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ne2</name>\n\
    <value>[O,S;X1]=*[$Ne1]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ne3</name>\n\
    <value>[nX2][n-]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ne4</name>\n\
    <value>[nX2][nX2][n-]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ne5</name>\n\
    <value>[nX2][nX2][nX2][n-]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ne6</name>\n\
    <value>[O,S;X1]=[$(ccc[A-]),$(CC=C[A-])]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Po0</name>\n\
    <value>[*;+]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Po1</name>\n\
    <value>[N,n;+;H,H2,H3]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Po2</name>\n\
    <value>[A;H,H2]-C=[$Po1]</value>\n\
    <vector_binding/>\n\
  </smarts>\n\
  <!-- final SMARTS used in pharmacophore features. Requires vector bindings. -->\n\
  <smarts>\n\
    <name>AllAcc</name>\n\
    <value>[$OAc,$NAc]</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>PosDon</name>\n\
    <value>[$Po1,$Po2]</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>NegAcc</name>\n\
    <value>[$Ne1,$Ne2,$Ne3,$Ne4,$Ne5,$Ne6;!$Ylid]</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>NeuDon</name>\n\
    <value>[$AllDon;!$PosDon]</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>NeuAcc</name>\n\
    <value>[$AllAcc;!$NegAcc]</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>AromOne</name>\n\
    <value>[aD2][aD3][aD2]</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>AromTwo</name>\n\
    <value>[aD2][aD2][aD2]</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>AromFive</name>\n\
    <value>a1aaaa1</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>AromSix</name>\n\
    <value>a1aaaaa1</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Ani</name>\n\
    <value>[$Ne1]</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <smarts>\n\
    <name>Cat</name>\n\
    <value>[$Po0]</value>\n\
    <final_definition/>\n\
  </smarts>\n\
  <!-- And finally, the features themselves. -->\n\
  <feature>\n\
    <name>Don</name>\n\
    <smarts>NeuDon</smarts>\n\
    <direction>h_vector</direction>\n\
  </feature>\n\
  <feature>\n\
    <name>Acc</name>\n\
    <smarts>NeuAcc</smarts>\n\
    <direction>lp_vector</direction>\n\
  </feature>\n\
  <feature>\n\
    <name>Ani</name>\n\
    <smarts>Ani</smarts>\n\
    <direction>lp_vector</direction>\n\
  </feature>\n\
  <feature>\n\
    <name>Cat</name>\n\
    <smarts>Cat</smarts>\n\
    <direction>h_vector</direction>\n\
  </feature>\n\
  <feature>\n\
    <name>Aro</name>\n\
    <smarts>AromFive</smarts>\n\
    <smarts>AromSix</smarts>\n\
    <direction>normal</direction>\n\
  </feature>\n\
  <feature>\n\
    <name>Hyd</name>\n\
    <algorithmic>aliphatic_only</algorithmic>\n\
    <direction>false</direction>\n\
  </feature>\n\
  <feature>\n\
    <name>AllHyd</name>\n\
    <algorithmic>aliphatic_only</algorithmic>\n\
    <smarts>AromFive</smarts>\n\
    <smarts>AromSix</smarts>\n\
    <direction>normal</direction>\n\
  </feature>\n\
</features>\n";
#ifdef NOTYET
#endif

} // end of namespace DACLIB

#endif // MOEPOINTSDEFS_H
