(* ::Package:: *)

Quit[]


$LoadTARCER = True;
$LoadFeynArts = True;
<< FeynCalc/FeynCalc.m;
AppendTo[$Path, 
  "/Users/jamesmckay/Documents/Programs/Mass_builder/src/"];
<< MassBuilder.m;
SetOptions[DiracSlash, Dimension -> D, 
 FeynCalcInternal -> True]; SetOptions[DiracTrace, 
 DiracTraceEvaluate -> True]; null = 0;


(* Determine one-loop counter-term couplings *)

path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa = 1/(16 \[Pi]^2);
(*sw=Sin[\[Theta]];
cw=Cos[\[Theta]];*)
(*F11*)

F11SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_1_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_2_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_3_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

F11SEdiv = makeFiniteAmplitude[F11SE, -1, D];

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_1_1c.mx"}]];

F11ct = makeFiniteAmplitude[SelfEnergyFinite*kappa, -1, D] /. 
    MassBuilderCTM1 -> 0 /. MassBuilderCTZ1 -> 0;

eq1 = FullSimplify[Coefficient[F11ct + F11SEdiv, p]];
eq2 = FullSimplify[Coefficient[F11ct + F11SEdiv, p, 0]];
solF11 = Solve[{eq1 == 0, eq2 == 0}, {d2Z, d2m}]

(*F12*)

F12SE = 0;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F12_g1_1_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F12_g1_2_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F12_g1_3_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SEdiv = makeFiniteAmplitude[F12SE, -1, D];
Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F12_g1_1_1c.mx"}]];
F12ct = makeFiniteAmplitude[SelfEnergyFinite*kappa, -1, D] /. 
    MassBuilderCTM1 -> 0 /. MassBuilderCTZ1 -> 0;
eq1 = FullSimplify[Coefficient[F12ct + F12SEdiv, p]];
eq2 = FullSimplify[Coefficient[F12ct + F12SEdiv, p, 0]];
solF12 = Solve[{eq1 == 0, eq2 == 0}, {d1Z, d1m}]
Set @@@ solF11[[1]];
Set @@@ solF12[[1]];


F11SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_1_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_2_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_3_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_1_1c.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;
F11SE = F11SE /. MassBuilderCTZ1 -> 0 /. MassBuilderCTZ2 -> 0 /. 
     MassBuilderCTM1 -> 0 /. MassBuilderCTM2 -> 0 /. 
   D -> 4 - 2 MassBuilderEpsilon;

C0 = Coefficient[FullSimplify[F11SE], MassBuilderEpsilon, 0];
C1 = Coefficient[FullSimplify[F11SE], MassBuilderEpsilon, 1];
C2 = Coefficient[FullSimplify[F11SE], MassBuilderEpsilon, 2];
SEtest = C0 + C1*MassBuilderEpsilon + C2*MassBuilderEpsilon^2;
FullSimplify[SEtest - F11SE];
SEN = FullSimplify[SEtest] /. MassBuilderEpsilon -> 0 /. Pair[Momentum[p], Momentum[p]] -> p^2;
SEN1 = SEN /. TBI[4, p^2, {{1, mw}, {1, MChi}}] -> Bwc[s]\
/. TAI[4, 0, {{1, mw}}] -> Aw /. TAI[4, 0, {{1, MChi}}] -> Ac\
/. TAI[4, 0, {{1, ma}}] -> Aa /. TAI[4, 0, {{1, mz}}] -> Az;

SEN1=SEN1/.p^2->s/.p->Sqrt[s];
DN = D[SEN1, s];
totalN1=2*MChi*FullSimplify[DN*SEN1];


(* express in terms of B1 *)
SEN2 = SEN /. TBI[4, p^2, {{1, mw}, {1, MChi}}] -> Bwc[s]\
/. TAI[4, 0, {{1, mw}}] -> Aw /. TAI[4, 0, {{1, MChi}}] -> Ac\
/. TAI[4, 0, {{1, ma}}] -> Aa /. TAI[4, 0, {{1, mz}}] -> Az;
SEN2 = FullSimplify[SEN2 /. (Aw) -> 2 p^2 B1cw[s] -(p^2 + MChi^2 - mw^2)*Bwc[s] + Ac];
Coefficient[SEN2 /. Pair[Momentum[p, D], Momentum[p, D]] -> p^2, p, -1];
SigmaKN = Coefficient[SEN2 /. Pair[Momentum[p, D], Momentum[p, D]] -> p^2, p,1];
SigmaMN = Coefficient[SEN2 /. Pair[Momentum[p, D], Momentum[p, D]] -> p^2, p,0];
DSigmaMN = D[SigmaMN, s];
DSigmaKN = D[SigmaKN, s];
totalN2 = FullSimplify[(SigmaMN + MChi*SigmaKN)*(SigmaKN + 2*MChi*DSigmaMN + 2*MChi^2*DSigmaKN)];

deltaM = totalN1-totalN2;

deltaM=Simplify[deltaM]/.MassBuilderEpsilon->0;
deltaMCForm=Simplify[deltaM\
/.Derivative[1][B1cw][s]->(1/(2*s))*(-2*B1cw[s]+Bcw[s]+(s+MChi^2-mw^2)*dBwc)\
/.Derivative[1][B1cz][s]->(1/(2*s))*(-2*B1cz[s]+Bcz[s]+(s+MChi^2-mz^2)*dBzc)\
/.Derivative[1][B1ca][s]->(1/(2*s))*(-2*B1ca[s]+Bca[s]+(s+MChi^2-ma^2)*dBac)\
/.B1cz[s]->(1/(2*s))*(Az-Ac+(s+MChi^2-mz^2)*Bzc[s])\
/.B1ca[s]->(1/(2*s))*(Aa-Ac+(s+MChi^2-ma^2)*Bac[s])\
/.B1cw[s]->(1/(2*s))*(Aw-Ac+(s+MChi^2-mw^2)*Bwc[s])\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bzc/.Bca[s]->Bac\
/.Bac[s]->Bac/.Bwc[s]->Bwc\
/.Bcw[s]->Bwc/.Bcz[s]->Bzc];
Simplify[deltaMCForm/.Sqrt[s]->MChi/.s->MChi^2]


(* F12 *)
F12SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_1_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_2_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_3_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_1_1c.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SE = F12SE/.MassBuilderCTZ1->0/.MassBuilderCTZ2->0/.MassBuilderCTM1->0/.MassBuilderCTM2->0/.D->4-2MassBuilderEpsilon;
C0=Coefficient[Simplify[F12SE],MassBuilderEpsilon,0];
C1=Coefficient[Simplify[F12SE],MassBuilderEpsilon,1];
C2=Coefficient[Simplify[F12SE],MassBuilderEpsilon,2];
SEtest2 = C0+C1*MassBuilderEpsilon+C2*MassBuilderEpsilon^2;
Simplify[SEtest2-F12SE]
SEC=Simplify[SEtest2]/.MassBuilderEpsilon->0/.Pair[Momentum[p],Momentum[p]]->p^2;

SEC1=SEC/.TBI[4,p^2,{{1,mz},{1,MChi}}]->Bzc[s]\
/.TBI[4,p^2,{{1,MChi},{1,ma}}]->Bac[s]\
/.TBI[4,p^2,{{1,mw},{1,MChi}}]->Bwc[s]\
/.TAI[4,0,{{1,mw}}]->Aw\
/.TAI[4,0,{{1,MChi}}]->Ac\
/.TAI[4,0,{{1,ma}}]->Aa\
/.TAI[4,0,{{1,mz}}]->Az;

SEC1=SEC1/.p^2->s/.p->Sqrt[s];
DC = D[SEC1, s];
totalC1=2*MChi*Simplify[DC*SEC1];


SEC2=SEC/.TBI[4,p^2,{{1,mz},{1,MChi}}]->Bzc[s]\
/.TBI[4,p^2,{{1,MChi},{1,ma}}]->Bac[s]\
/.TBI[4,p^2,{{1,mw},{1,MChi}}]->Bwc[s]\
/.TAI[4,0,{{1,mw}}]->Aw\
/.TAI[4,0,{{1,MChi}}]->Ac\
/.TAI[4,0,{{1,ma}}]->Aa\
/.TAI[4,0,{{1,mz}}]->Az;
SEC2=Simplify[SEC2\
/. (Aw) -> 2 p^2 B1cw[s] -(p^2 + MChi^2 - mw^2)*Bwc[s] + Ac\
/. (Az) -> 2 p^2 B1cz[s] -(p^2 + MChi^2 - mz^2)*Bzc[s] + Ac\
/. (Aa) -> 2 p^2 B1ca[s] -(p^2 + MChi^2 - ma^2)*Bac[s] + Ac];

Coefficient[SEC2 /. Pair[Momentum[p, D], Momentum[p, D]] -> p^2, p, -1]
SigmaKC = Coefficient[SEC2 /. Pair[Momentum[p, D], Momentum[p, D]] -> p^2, p,1]
SigmaMC = Coefficient[SEC2 /. Pair[Momentum[p, D], Momentum[p, D]] -> p^2, p,0]
DSigmaMC = D[SigmaMC, s];
DSigmaKC = D[SigmaKC, s];
totalC2 = FullSimplify[(SigmaMC + MChi*SigmaKC)*(SigmaKC + 2*MChi*DSigmaMC + 2*MChi^2*DSigmaKC)];


deltaM = totalC1-totalC2;

deltaM=Simplify[deltaM]/.MassBuilderEpsilon->0;
deltaMCForm=Simplify[deltaM\
/.Derivative[1][B1cw][s]->(1/(2*s))*(-2*B1cw[s]+Bcw[s]+(s+MChi^2-mw^2)*dBwc)\
/.Derivative[1][B1cz][s]->(1/(2*s))*(-2*B1cz[s]+Bcz[s]+(s+MChi^2-mz^2)*dBzc)\
/.Derivative[1][B1ca][s]->(1/(2*s))*(-2*B1ca[s]+Bca[s]+(s+MChi^2-ma^2)*dBac)\
/.B1cz[s]->(1/(2*s))*(Az-Ac+(s+MChi^2-mz^2)*Bzc[s])\
/.B1ca[s]->(1/(2*s))*(Aa-Ac+(s+MChi^2-ma^2)*Bac[s])\
/.B1cw[s]->(1/(2*s))*(Aw-Ac+(s+MChi^2-mw^2)*Bwc[s])\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bzc/.Bca[s]->Bac\
/.Bac[s]->Bac/.Bwc[s]->Bwc\
/.Bcw[s]->Bwc/.Bcz[s]->Bzc]/.Sqrt[s]->MChi/.s->MChi^2;
del=Simplify[deltaMCForm];
FullSimplify[del,{MChi>0}]


finaldeltaM=Simplify[ (totalC1-totalN1) - (totalC2-totalN2)]/.MassBuilderEpsilon->0;
finaldeltaM2=Simplify[finaldeltaM\
/.Derivative[1][B1cw][s]->(1/(2*s))*(-2*B1cw[s]+Bcw[s]+(s+MChi^2-mw^2)*dBwc)\
/.Derivative[1][B1cz][s]->(1/(2*s))*(-2*B1cz[s]+Bcz[s]+(s+MChi^2-mz^2)*dBzc)\
/.Derivative[1][B1ca][s]->(1/(2*s))*(-2*B1ca[s]+Bca[s]+(s+MChi^2-ma^2)*dBac)\
/.B1cz[s]->(1/(2*s))*(Az-Ac+(s+MChi^2-mz^2)*Bzc[s])\
/.B1ca[s]->(1/(2*s))*(Aa-Ac+(s+MChi^2-ma^2)*Bac[s])\
/.B1cw[s]->(1/(2*s))*(Aw-Ac+(s+MChi^2-mw^2)*Bwc[s])\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bcz/.Bca[s]->Bca\
/.Bac[s]->Bca/.Bwc[s]->Bcw\
/.Bcw[s]->Bcw/.Bcz[s]->Bcz]/.Sqrt[s]->MChi/.s->MChi^2;

Simplify[finaldeltaM2,{MChi>0}]


finaldeltaM=Simplify[ (totalC1-totalN1)]/.MassBuilderEpsilon->0;
finaldeltaM2=Simplify[finaldeltaM\
/.Derivative[1][B1cw][s]->(1/(2*s))*(-2*B1cw[s]+Bcw[s]+(s+MChi^2-mw^2)*dBwc)\
/.Derivative[1][B1cz][s]->(1/(2*s))*(-2*B1cz[s]+Bcz[s]+(s+MChi^2-mz^2)*dBzc)\
/.Derivative[1][B1ca][s]->(1/(2*s))*(-2*B1ca[s]+Bca[s]+(s+MChi^2-ma^2)*dBac)\
/.B1cz[s]->(1/(2*s))*(Az-Ac+(s+MChi^2-mz^2)*Bzc[s])\
/.B1ca[s]->(1/(2*s))*(Aa-Ac+(s+MChi^2-ma^2)*Bac[s])\
/.B1cw[s]->(1/(2*s))*(Aw-Ac+(s+MChi^2-mw^2)*Bwc[s])\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bcz/.Bca[s]->Bca\
/.Bac[s]->Bca/.Bwc[s]->Bcw\
/.Bcw[s]->Bcw/.Bcz[s]->Bcz]/.Sqrt[s]->MChi/.s->MChi^2;

CForm[Simplify[finaldeltaM2,{MChi>0}]]
