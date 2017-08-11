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
    "/models/MSSM/output/math_data_F11_g1_1_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/MSSM/output/math_data_F11_g1_2_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/MSSM/output/math_data_F11_g1_3_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

F11SEdiv = makeFiniteAmplitude[F11SE, -1, D];

Get[FileNameJoin[{path, 
    "/models/MSSM/output/math_data_F11_g1_1_1c.mx"}]];

F11ct = makeFiniteAmplitude[SelfEnergyFinite*kappa, -1, D] /. 
    MassBuilderCTM1 -> 0 /. MassBuilderCTZ1 -> 0;

eq1 = FullSimplify[Coefficient[F11ct + F11SEdiv, p]];
eq2 = FullSimplify[Coefficient[F11ct + F11SEdiv, p, 0]];
solF11 = Solve[{eq1 == 0, eq2 == 0}, {d2Z, d2m}]

(*F12*)

F12SE = 0;

Get[FileNameJoin[{path, 
    "/models/MSSM/output/math_data_F12_g1_1_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/MSSM/output/math_data_F12_g1_2_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/MSSM/output/math_data_F12_g1_3_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SEdiv = makeFiniteAmplitude[F12SE, -1, D];
Get[FileNameJoin[{path, 
    "/models/MSSM/output/math_data_F12_g1_1_1c.mx"}]];
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
    "/models/MSSM/output/math_data_F11_g1_1_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/MSSM/output/math_data_F11_g1_2_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/MSSM/output/math_data_F11_g1_3_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/MSSM/output/math_data_F11_g1_1_1c.mx"}]];
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
SEN1 = SEN /. TBI[4, p^2, {{1, mw}, {1, MChi}}] ->b I*Bwc[s]\
/. TAI[4, 0, {{1, mw}}] -> a*I*Aw \
/. TAI[4, 0, {{1, MChi}}] -> a*I*Ac\
/. TAI[4, 0, {{1, ma}}] -> a*I*Aa \
/. TAI[4, 0, {{1, mz}}] -> a*I*Az;

SEN1=SEN1/.p^2->s/.p->Sqrt[s];
DN = D[SEN1, s];
totalN1=2*MChi*FullSimplify[DN*SEN1];
(* F12 *)
F12SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_1_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_2_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_3_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_1_1c.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SE = F12SE/.MassBuilderCTZ1->0/.MassBuilderCTZ2->0/.MassBuilderCTM1->0/.MassBuilderCTM2->0/.D->4-2MassBuilderEpsilon;
C0=Coefficient[Simplify[F12SE],MassBuilderEpsilon,0];
C1=Coefficient[Simplify[F12SE],MassBuilderEpsilon,1];
C2=Coefficient[Simplify[F12SE],MassBuilderEpsilon,2];
SEtest2 = C0+C1*MassBuilderEpsilon+C2*MassBuilderEpsilon^2;
Simplify[SEtest2-F12SE]
SEC=Simplify[SEtest2]/.MassBuilderEpsilon->0/.Pair[Momentum[p],Momentum[p]]->p^2;

SEC1=SEC/.TBI[4,p^2,{{1,mz},{1,MChi}}]->b I*Bzc[s]\
/.TBI[4,p^2,{{1,MChi},{1,ma}}]->I b*Bac[s]\
/.TBI[4,p^2,{{1,mw},{1,MChi}}]->I b *Bwc[s]\
/.TAI[4,0,{{1,mw}}]->a I*Aw\
/.TAI[4,0,{{1,MChi}}]->a I*Ac\
/.TAI[4,0,{{1,ma}}]->a I*Aa\
/.TAI[4,0,{{1,mz}}]->a I*Az;

SEC1=SEC1/.p^2->s/.p->Sqrt[s];
DC = D[SEC1, s];
totalC1=2*MChi*Simplify[DC*SEC1];


a=1;
b=1;
deltaM= (totalN1-totalC1)/.MassBuilderEpsilon->0;
deltaM2=Simplify[deltaM\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bcz/.Bca[s]->Bca\
/.Bac[s]->Bca/.Bwc[s]->Bcw\
/.Bcw[s]->Bcw/.Bcz[s]->Bcz]/.Sqrt[s]->MChi/.s->MChi^2;

CForm[Simplify[deltaM2,{MChi>0}]]


(* Ibe et al. working *)
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];
(*Clear[cw,sw];*)
SigmaKN=-(e^2/(16*Pi^2*sw^2))*(4*B1cw[s]+2);
SigmaMN=-(e^2*MChi/(16*Pi^2*sw^2))*(8*Bcw[s]-4);
DSigmaMN=D[SigmaMN,s];
DSigmaKN=D[SigmaKN,s];
totalN =Simplify[(SigmaMN+MChi*SigmaKN)*(SigmaKN+2*MChi*DSigmaMN+2*MChi^2*DSigmaKN)];
SigmaKC = -(e^2/(8*Pi^2 * sw^2)) * ( sw^2* B1ca[s] + cw^2 *B1cz[s] + B1cw[s] + 1 );
SigmaMC = -(e^2 MChi /(4*Pi^2 * sw^2)) * ( sw^2* Bca[s] +  cw^2 *Bcz[s]+  Bcw[s] - 1 );
DSigmaMC = D[SigmaMC,s];
DSigmaKC = D[SigmaKC,s];
totalC=Simplify[(SigmaMC+MChi*SigmaKC)*(SigmaKC+2*MChi*DSigmaMC+2*MChi^2*DSigmaKC)];
deltaM=totalN-totalC/.MassBuilderEpsilon->0;
deltaMIbe=Simplify[deltaM\
/.Derivative[1][B1cw][s]-> (1/(2*s))*(-2*B1cw[s]-Bcw[s]-(s+MChi^2-mw^2)*dBwc)\
/.Derivative[1][B1cz][s]->(1/(2*s))*(-2*B1cz[s]-Bcz[s]-(s+MChi^2-mz^2)*dBzc)\
/.Derivative[1][B1ca][s]->(1/(2*s))*(-2*B1ca[s]-Bca[s]-(s+MChi^2-ma^2)*dBac)\
/.B1cz[s]->(1/(2*s))*(Ac-Az-(s+MChi^2-mz^2)*Bzc[s])\
/.B1ca[s]->(1/(2*s))*(Ac-Aa-(s+MChi^2-ma^2)*Bac[s])\
/.B1cw[s]->(1/(2*s))*(Ac-Aw-(s+MChi^2-mw^2)*Bwc[s])\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bcz/.Bca[s]->Bca\
/.Bac[s]->Bca/.Bwc[s]->Bcw\
/.Bcw[s]->Bcw/.Bcz[s]->Bcz]/.Sqrt[s]->MChi/.s->MChi^2;


a=1;
b=1;
Simplify[Coefficient[deltaMIbe-deltaM2,dBac],{MChi>0}]



FullSimplify[deltaMIbe-deltaM2,{MChi>0}]


CForm[deltaMIbe]


a=1;
b=1;
deltaM= (totalN1-totalC1)/.MassBuilderEpsilon->0;
totalN1CForm=Simplify[totalN1\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bcz/.Bca[s]->Bca\
/.Bac[s]->Bca/.Bwc[s]->Bcw\
/.Bcw[s]->Bcw/.Bcz[s]->Bcz]/.Sqrt[s]->MChi/.s->MChi^2;

totalC1CForm=Simplify[totalC1\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bcz/.Bca[s]->Bca\
/.Bac[s]->Bca/.Bwc[s]->Bcw\
/.Bcw[s]->Bcw/.Bcz[s]->Bcz]/.Sqrt[s]->MChi/.s->MChi^2;


CForm[Simplify[totalN1CForm,{MChi>0}]]


CForm[Simplify[totalC1CForm,{MChi>0}]]
