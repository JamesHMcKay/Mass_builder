(* ::Package:: *)

Quit[]


(* ::Section:: *)
(*Determine counter-term couplings for diagrams in group B*)


(* ::Subsection:: *)
(* Load required packages *)


$LoadTARCER = True;
$LoadFeynArts = True;
<< FeynCalc/FeynCalc.m;
AppendTo[$Path, "/Users/jamesmckay/Documents/Programs/Mass_builder/src/"];
<<MassBuilder.m;


SetOptions[DiracSlash,Dimension->D,FeynCalcInternal->True];SetOptions[DiracTrace,DiracTraceEvaluate->True];null=0;MassBuilderA[mass_,D_]:=TAI[D,0,{{1,mass}}];MassBuilderB[mass1_,mass2_,D_]:=TBI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2}}];MassBuilderJ[mass1_,mass2_,mass3_,D_]:=TJI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3}}];MassBuilderT[mass1_,mass2_,mass3_,D_]:=TJI[D,Pair[Momentum[p,D],Momentum[p,D]],{{2,mass1},{1,mass2},{1,mass3}}];MassBuilderK[mass1_,mass2_,mass3_,D_]:=TJI[D,0,{{1,mass1},{1,mass2},{1,mass3}}];MassBuilderV[mass1_,mass2_,mass3_,mass4_,D_]:=TVI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3},{1,mass4}}];MassBuilderF[mass1_,mass2_,mass3_,mass4_,mass5_,D_]:=TFI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3},{1,mass4},{1,mass5}}];FCGV["EL"]=e;FCGV["SW"]=sw;FCGV["CW"]=cw;FCGV["MW"]=mw;FCGV["MZ"]=mz;FCGV["ME"]=me;FCGV["MM"]=mm;FCGV["ML"]=ml;FCGV["MU"]=mu;FCGV["MT"]=mt;FCGV["MD"]=md;FCGV["MS"]=ms;FCGV["MB"]=mb;ZNeu[1,1]=0;ZNeu[1,2]=1;ZNeu[1,3]=0;ZNeu[1,4]=0;ZNeu[2,1]=0;ZNeu[2,2]=0;ZNeu[2,3]=0;ZNeu[2,4]=0;ZNeu[3,1]=0;ZNeu[3,2]=0;ZNeu[3,3]=0;ZNeu[3,4]=0;ZNeu[4,1]=0;ZNeu[4,2]=0;ZNeu[4,3]=0;ZNeu[4,4]=0;UCha[1,1]=1;UCha[1,2]=0;UCha[2,1]=0;UCha[2,2]=0;VCha[1,1]=1;VCha[1,2]=0;VCha[2,1]=0;VCha[2,2]=0;CKM=IndexDelta;MCha[1]=MChi;MCha[2]=MChi;MNeu[1]=MChi;MNeu[2]=MChi;Mh0=mh;SBA=1;


(* ::Subsection:: *)
(* Get the amplitudes that have been computed previously*)
(* Change the path below to the location of your Mass_Builder/models directory *)


path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa=1/(16\[Pi]^2);
STW=Sin[\[Theta]];
CTW=Cos[\[Theta]];

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


SPD[p,p]=MChi^2;Pair[Momentum[p],Momentum[p]]=MChi^2;
(* Particle 1 *)
SelfEnergyTotal = 0;


Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_131_2.mx"}]];
C131 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_132_2.mx"}]];
C132 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_133_2.mx"}]];
C133 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_134_2.mx"}]];
C134 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_135_2.mx"}]];
C135 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_136_2.mx"}]];
C136 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_137_2.mx"}]];
C137 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_138_2.mx"}]];
C138 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_139_2.mx"}]];
C139 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_140_2.mx"}]];
C140 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_141_2.mx"}]];
C141 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_142_2.mx"}]];
C142 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_143_2.mx"}]];
C143 = SelfEnergyFinite*kappa^2;
Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_1_2c.mx"}]];
C1c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_2_2c.mx"}]];
C2c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_3_2c.mx"}]];
C3c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_4_2c.mx"}]];
C4c = SelfEnergyFinite*kappa^2;

Csum1=makeFiniteAmplitude[C131+C132+C133+C134+C135+C136+C137+C138+C139+C140+C141+C142+C143,-1,D];
Csum2=makeFiniteAmplitude[C131+C132+C133+C134+C135+C136+C137+C138+C139+C140+C141+C142+C143,-2,D];

Cct1 = makeFiniteAmplitude[C1c+C2c+C3c+C4c,-1,D];
Cct2 = makeFiniteAmplitude[C1c+C2c+C3c+C4c,-2,D];



(* Particle 2 *)

SelfEnergyTotal = 0;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_99_2.mx"}]];
CC99 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_100_2.mx"}]];
CC100 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_101_2.mx"}]];
CC101 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_102_2.mx"}]];
CC102 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_103_2.mx"}]];
CC103 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_104_2.mx"}]];
CC104 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_105_2.mx"}]];
CC105 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_106_2.mx"}]];
CC106 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_107_2.mx"}]];
CC107 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_108_2.mx"}]];
CC108 = SelfEnergyFinite*kappa^2;


Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_1_2c.mx"}]];
CC1c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_2_2c.mx"}]];
CC2c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_3_2c.mx"}]];
CC3c = SelfEnergyFinite*kappa^2;

CCsum1=makeFiniteAmplitude[CC99+CC100+CC101+CC102+CC103+CC104+CC105+CC106+CC107+CC108,-1,D];
CCsum2=makeFiniteAmplitude[CC99+CC100+CC101+CC102+CC103+CC104+CC105+CC106+CC107+CC108,-2,D];

CCct1 = makeFiniteAmplitude[CC1c+CC2c+CC3c,-1,D];
CCct2 = makeFiniteAmplitude[CC1c+CC2c+CC3c,-2,D];


(* Particle 3 *)

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_61_2.mx"}]];
N61 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_62_2.mx"}]];
N62 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_63_2.mx"}]];
N63 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_64_2.mx"}]];
N64 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_65_2.mx"}]];
N65 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_66_2.mx"}]];
N66 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_67_2.mx"}]];
N67 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_68_2.mx"}]];
N68 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_1_2c.mx"}]];
N1c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_2_2c.mx"}]];
N2c = SelfEnergyFinite*kappa^2;

Nsum1=makeFiniteAmplitude[N61+N62+N63+N64+N65+N66+N67+N68,-1,D];
Nsum2=makeFiniteAmplitude[N61+N62+N63+N64+N65+N66+N67+N68,-2,D];

Nct1 = makeFiniteAmplitude[N1c+N2c,-1,D];
Nct2 = makeFiniteAmplitude[N1c+N2c,-2,D];




{Csum1,Csum2,Nsum1,Nsum2,Cct1,Cct2,Nct1,Nct2,CCsum1,CCsum2,CCct1,CCct2}={Csum1,Csum2,Nsum1,Nsum2,Cct1,Cct2,Nct1,Nct2,CCsum1,CCsum2,CCct1,CCct2}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
/. MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]]/. MassBuilderP -> Momentum[p]/. MassBuilderEpsilon -> \[Epsilon]/. MassBuilderZeta -> \[Zeta]/. MassBuilderAe[0] -> 0\
/. D -> 4 - 2 \[Epsilon]/. Pair[p, p] -> p^2/. D -> 4 - 2 \[Epsilon]/. Momentum[p] -> p\
/. MassBuilderB[MChi, mz] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[mz, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[MChi, mw] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[mw, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[ma, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, ma}}]\
/. MassBuilderB[MChi, ma] -> TBI[4, MChi^2, {{1, MChi}, {1, ma}}]\
/.p->MChi/.STW->sw/.CTW->cw;
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];


FullSimplify[Nsum2+Nct2-(Csum2+Cct2)/.mz->mw/cw]
FullSimplify[Nsum2+Nct2-(CCsum2+CCct2)/.mz->mw/cw]


FullSimplify[Nsum2+Nct2/.mz->mw/cw]
FullSimplify[Csum2+Cct2/.mz->mw/cw]
FullSimplify[CCsum2+CCct2/.mz->mw/cw]


FullSimplify[Nsum1+Nct1-(Csum1+Cct1)/.mz->mw/cw]
FullSimplify[Nsum1+Nct1-(CCsum1+CCct1)/.mz->mw/cw]


CForm[FullSimplify[d1Z]]


CForm[FullSimplify[d1m]]


CForm[FullSimplify[d2Z]]


CForm[FullSimplify[d2m]]


FullSimplify[d5m*kappa]
FullSimplify[d5m*kappa]
FullSimplify[d6Z*kappa]
FullSimplify[d6Z*kappa]
FullSimplify[d7Z*kappa]
FullSimplify[d7Z*kappa]



