#
#  Proposed SMARTS for Plurality
#  
#  It is intended that these SMARTS be used with a suitable
#  ionisation and tautomer model.
#
#  Hydrogen bond donors and acceptors are defined first
#  without reference to charge.  Charged donors and 
#  acceptors are defined for use in their own right and to
#  specify neutral acceptors and donors. Ylids are defined
#  as a special case because CORINA converts pentavalent
#  representation of hypervalent nitrogen to ylid in the
#  generation of 3D coordinates. 
#
#----------------------------------------------------------
#  Define ylids beacuse CORINA will convert pentavalent
#  nitrogen form of hypervalent species to ylid form
Ylid    [A-]~[N,n;+;H0]                                0  0
#----------------------------------------------------------
# Define donors (neutral & cationic) 
AllDon  [O,N,n;!H0]                                    1  1
#----------------------------------------------------------
#  Define acceptors (neutral & anionic)
Ac1     [O;X1]                                         1  0
Ac2     [OH]                                           1  0
Ac3     O([CX4])[CX4]                                  1  0
OAc     [$Ac1,$Ac2,$Ac3]                               1  0
Ac4     [nX2]                                          1  0
Ac5     [NX2]=C                                        1  0
Ac6     [N;X3;H2][CX4]                                 1  0
Ac7     [N;X3;H]([CX4])[CX4]                           1  0
Ac8     [N;X3]([CX4])([CX4])[CX4]                      1  0
Ac9     [N;X1]                                         1  0
# Adding in SMARTS to take out none acceptors in oxadiazoles
NonAccN  [$(ncC=O)]                                    1  0
NAc     [$Ac4,$Ac5,$Ac6,$Ac7,$Ac8,$Ac9;!$NonAccN]                1  0
AllAcc  [$OAc,$NAc]                                    1  1
#----------------------------------------------------------
# Define catonic donors
Po1     [N,n;+;H,H2,H3]                                1  0
Po2     [A;H,H2]-C=[$Po1]                              1  0
PosDon  [$Po1,$Po2]                                    1  1
#----------------------------------------------------------
# Define anionic acceptors
Ne1     [*;-]                                          1  0
Ne2     [O,S;X1]=*[$Ne1]                               1  0
Ne3     [nX2][n-]                                      1  0
Ne4     [nX2][nX2][n-]                                 1  0
Ne5     [nX2][nX2][nX2][n-]                            1  0
Ne6     [O,S;X1]=[$(ccc[A-]),$(CC=C[A-])]              1  0
NegAcc  [$Ne1,$Ne2,$Ne3,$Ne4,$Ne5,$Ne6;!$Ylid]         1  1
#----------------------------------------------------------
# Define neutral donors & acceptors by eliminating 
# charged species from all donors & acceptors
NeuDon  [$AllDon;!$PosDon]                             1  1
NeuAcc  [$AllAcc;!$NegAcc]                             1  1
#----------------------------------------------------------
AromOne  [aD2][aD3][aD2]                                      1  1
AromTwo   [aD2][aD2][aD2]                                     1  1
AromFive a1aaaa1					1 1
AromSix  a1aaaaa1					1 1

