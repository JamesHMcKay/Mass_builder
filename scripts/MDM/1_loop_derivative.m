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
kappa=1/(16\[Pi]^2);
(*STW=Sin[\[Theta]];
CTW=Cos[\[Theta]];*)

(* F5 *)

F5SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_1_1.mx"}]];
F5SE = F5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_2_1.mx"}]];
F5SE = F5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_3_1.mx"}]];
F5SE = F5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_4_1.mx"}]];
F5SE = F5SE + SelfEnergyFinite*kappa;

F5SEdiv = makeFiniteAmplitude[F5SE,-1,D];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_1_1c.mx"}]];

F5ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;

eq1 = FullSimplify[Coefficient[F5ct+F5SEdiv,p]];
eq2 = FullSimplify[Coefficient[F5ct+F5SEdiv,p,0]];
solF5 = Solve[{eq1==0,eq2==0},{d5Z, d5m}]

(* F6 *)

F6SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_1_1.mx"}]];
F6SE = F6SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_2_1.mx"}]];
F6SE = F6SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_3_1.mx"}]];
F6SE = F6SE + SelfEnergyFinite*kappa;
F6SEdiv = makeFiniteAmplitude[F6SE,-1,D];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_1_1c.mx"}]];

F6ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;

eq1 = FullSimplify[Coefficient[F6ct+F6SEdiv,p]];
eq2 = FullSimplify[Coefficient[F6ct+F6SEdiv,p,0]];
solF6 = Solve[{eq1==0,eq2==0},{d6Z, d6m}]


(* F7 *)

F7SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_1_1.mx"}]];
F7SE = F7SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_2_1.mx"}]];
F7SE = F7SE + SelfEnergyFinite*kappa;
F7SEdiv = makeFiniteAmplitude[F7SE,-1,D];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_1_1c.mx"}]];

F7ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;

eq1 = FullSimplify[Coefficient[F7ct+F7SEdiv,p]];
eq2 = FullSimplify[Coefficient[F7ct+F7SEdiv,p,0]];
solF7 = Solve[{eq1==0,eq2==0},{d7Z, d7m}]
Set @@@ solF5[[1]];
Set @@@ solF6[[1]];
Set @@@ solF7[[1]];


F11SE = 0;
ClearScalarProducts[];

(* F5 *)

F5SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_1_1.mx"}]];
F5SE = F5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_2_1.mx"}]];
F5SE = F5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_3_1.mx"}]];
F5SE = F5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_4_1.mx"}]];
F5SE = F5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_1_1c.mx"}]];
F5SE = F5SE + SelfEnergyFinite*kappa;

F5SE = F5SE /. MassBuilderCTZ1 -> 0 /. MassBuilderCTZ2 -> 0 /. 
     MassBuilderCTM1 -> 0 /. MassBuilderCTM2 -> 0 /. 
   D -> 4 - 2 MassBuilderEpsilon;

C0 = Coefficient[FullSimplify[F5SE], MassBuilderEpsilon, 0];
C1 = Coefficient[FullSimplify[F5SE], MassBuilderEpsilon, 1];
C2 = Coefficient[FullSimplify[F5SE], MassBuilderEpsilon, 2];
SEtest = C0 + C1*MassBuilderEpsilon + C2*MassBuilderEpsilon^2;
Simplify[SEtest - F5SE]


SEC = FullSimplify[SEtest] /. MassBuilderEpsilon -> 0 /. Pair[Momentum[p], Momentum[p]] -> p^2;
SEC1 = SEC /. TBI[4, p^2, {{1, mw}, {1, MChi}}] ->b I*Bwc[s]\
/.TBI[4,p^2,{{1,MChi},{1,ma}}]->I b*Bac[s]\
/.TBI[4,p^2,{{1,mw},{1,MChi}}]->I b *Bwc[s]\
/.TBI[4,p^2,{{1,mz},{1,MChi}}]->I b *Bzc[s]\
/. TAI[4, 0, {{1, mw}}] -> a*I*Aw \
/. TAI[4, 0, {{1, MChi}}] -> a*I*Ac\
/. TAI[4, 0, {{1, ma}}] -> a*I*Aa \
/. TAI[4, 0, {{1, mz}}] -> a*I*Az


Simplify[Coefficient[F5SE,MassBuilderEpsilon,-1]]



SEC1=SEC1/.p^2->s/.p->Sqrt[s];
DC = D[SEC1, s];
totalC1=2*MChi*FullSimplify[DC*SEC1];



(* F6 *)

F6SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_1_1.mx"}]];
F6SE = F6SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_2_1.mx"}]];
F6SE = F6SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_3_1.mx"}]];
F6SE = F6SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_1_1c.mx"}]];
F6SE = F6SE + SelfEnergyFinite*kappa;


F6SE = F6SE/.MassBuilderCTZ1->0/.MassBuilderCTZ2->0/.MassBuilderCTM1->0/.MassBuilderCTM2->0/.D->4-2MassBuilderEpsilon;
C0=Coefficient[Simplify[F6SE],MassBuilderEpsilon,0];
C1=Coefficient[Simplify[F6SE],MassBuilderEpsilon,1];
C2=Coefficient[Simplify[F6SE],MassBuilderEpsilon,2];
SEtest2 = C0+C1*MassBuilderEpsilon+C2*MassBuilderEpsilon^2;
Simplify[SEtest2-F6SE]
SECC=Simplify[SEtest2]/.MassBuilderEpsilon->0/.Pair[Momentum[p],Momentum[p]]->p^2;

SECC1=SECC/.TBI[4,p^2,{{1,mz},{1,MChi}}]->b I*Bzc[s]\
/.TBI[4,p^2,{{1,MChi},{1,ma}}]->I b*Bac[s]\
/.TBI[4,p^2,{{1,mw},{1,MChi}}]->I b *Bwc[s]\
/.TAI[4,0,{{1,mw}}]->a I*Aw\
/.TAI[4,0,{{1,MChi}}]->a I*Ac\
/.TAI[4,0,{{1,ma}}]->a I*Aa\
/.TAI[4,0,{{1,mz}}]->a I*Az;

SECC1=SECC1/.p^2->s/.p->Sqrt[s];
DCC = D[SECC1, s];
totalCC1=2*MChi*Simplify[DCC*SECC1];




Simplify[Coefficient[F6SE,MassBuilderEpsilon,-1]]


(* F7 *)

F7SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_1_1.mx"}]];
F7SE = F7SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_2_1.mx"}]];
F7SE = F7SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_1_1c.mx"}]];
F7SE = F7SE + SelfEnergyFinite*kappa;


F7SE = F7SE/.MassBuilderCTZ1->0/.MassBuilderCTZ2->0/.MassBuilderCTM1->0/.MassBuilderCTM2->0/.D->4-2MassBuilderEpsilon;
C0=Coefficient[Simplify[F7SE],MassBuilderEpsilon,0];
C1=Coefficient[Simplify[F7SE],MassBuilderEpsilon,1];
C2=Coefficient[Simplify[F7SE],MassBuilderEpsilon,2];
SEtest3 = C0+C1*MassBuilderEpsilon+C2*MassBuilderEpsilon^2;
Simplify[SEtest3-F7SE]
SEN=Simplify[SEtest3]/.MassBuilderEpsilon->0/.Pair[Momentum[p],Momentum[p]]->p^2;

SEN1=SEN/.TBI[4,p^2,{{1,mz},{1,MChi}}]->b I*Bzc[s]\
/.TBI[4,p^2,{{1,MChi},{1,ma}}]->I b*Bac[s]\
/.TBI[4,p^2,{{1,mw},{1,MChi}}]->I b *Bwc[s]\
/.TAI[4,0,{{1,mw}}]->a I*Aw\
/.TAI[4,0,{{1,MChi}}]->a I*Ac\
/.TAI[4,0,{{1,ma}}]->a I*Aa\
/.TAI[4,0,{{1,mz}}]->a I*Az;

SEN1=SEN1/.p^2->s/.p->Sqrt[s];
DN = D[SEN1, s];
totalN1=2*MChi*Simplify[DN*SEN1];


Simplify[Coefficient[F7SE,MassBuilderEpsilon,-1]]


a=1;
b=1;
deltaN1= (totalN1)/.MassBuilderEpsilon->0;
deltaN2=Simplify[deltaN1\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bcz/.Bca[s]->Bca\
/.Bac[s]->Bca/.Bwc[s]->Bcw\
/.Bcw[s]->Bcw/.Bcz[s]->Bcz]/.Sqrt[s]->MChi/.s->MChi^2;

CForm[Simplify[deltaN2,{MChi>0}]]


a=1;
b=1;
deltaC1= (totalC1)/.MassBuilderEpsilon->0;
deltaC2=Simplify[deltaC1\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bcz/.Bca[s]->Bca\
/.Bac[s]->Bca/.Bwc[s]->Bcw\
/.Bcw[s]->Bcw/.Bcz[s]->Bcz]/.Sqrt[s]->MChi/.s->MChi^2;
CForm[Simplify[deltaC2,{MChi>0}]]


a=1;
b=1;
deltaCC1= (totalCC1)/.MassBuilderEpsilon->0;
deltaCC2=Simplify[deltaCC1\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bcz/.Bca[s]->Bca\
/.Bac[s]->Bca/.Bwc[s]->Bcw\
/.Bcw[s]->Bcw/.Bcz[s]->Bcz]/.Sqrt[s]->MChi/.s->MChi^2;

CForm[Simplify[deltaCC2,{MChi>0}]]
