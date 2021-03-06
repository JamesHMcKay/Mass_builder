(* ::Package:: *)

Quit[]


(* ::Section:: *)
(*Determine counter-term couplings for diagrams in group C*)


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
(*sw=Sin[\[Theta]];
cw=Cos[\[Theta]];*)
V3SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_7_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_8_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_9_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_10_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_11_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_12_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

V3SEdiv = makeFiniteAmplitude[V3SE,-1,D];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_1_1c.mx"}]];
V3ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V3ct+V3SEdiv,Pair[Momentum[p],Momentum[p]]]];
eq2 = FullSimplify[Coefficient[V3ct+V3SEdiv,Pair[Momentum[p],Momentum[p]],0]];
solV3 = Solve[{eq1==0,eq2==0},{dZW1, dMWsq1}]
(* Z boson *) 

V2SE = 0;
ClearScalarProducts[];


Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_5_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_6_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_7_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_8_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_9_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_10_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_11_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_12_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_13_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_14_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_15_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_16_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;


V2SEdiv = makeFiniteAmplitude[V2SE,-1,D];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_1_1c.mx"}]];
V2ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V2ct+V2SEdiv,Pair[Momentum[p],Momentum[p]]]];
eq2 = FullSimplify[Coefficient[V2ct+V2SEdiv,Pair[Momentum[p],Momentum[p]],0]];
solV2 = Solve[{eq1==0,eq2==0},{dMZsq1,dZZZ1}]

(* photon *) 
V1SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_3_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_4_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_5_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_6_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_7_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_8_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_9_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_10_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_11_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;


V1SEdiv = makeFiniteAmplitude[V1SE,-1,D];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_1_1c.mx"}]];
V1ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V1ct+V1SEdiv,Pair[Momentum[p],Momentum[p]]]];
(*solV1 = Solve[{eq1\[Equal]0,eq2\[Equal]0},{dZAA1,dMAsq1}]*)
solV1 = Solve[{eq1==0},{dZAA1}]

(* photon Z mixing *) 

V12SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_3_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_4_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_5_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_6_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_7_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_8_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_9_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_10_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_11_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

V12SEdiv = makeFiniteAmplitude[V12SE,-1,D];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_1_1c.mx"}]];
V12ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V12ct+V12SEdiv,Pair[Momentum[p],Momentum[p]]]];
eq2 = FullSimplify[Coefficient[V12ct+V12SEdiv,Pair[Momentum[p],Momentum[p]],0]];
solV12 = Solve[{eq1==0,eq2==0},{dZZA1,dZAZ1}]



Set @@@ solV1[[1]];
Set @@@ solV2[[1]];
Set @@@ solV3[[1]];
Set @@@ solV12[[1]];


SPD[p,p]=MChi^2;Pair[Momentum[p],Momentum[p]]=MChi^2;
(* Particle 1 *)

SelfEnergyTotal = 0;


Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_19_2.mx"}]];
N19 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_20_2.mx"}]];
N20 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_21_2.mx"}]];
N21 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_22_2.mx"}]];
N22 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_23_2.mx"}]];
N23 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_24_2.mx"}]];
N24 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_25_2.mx"}]];
N25 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_26_2.mx"}]];
N26 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_27_2.mx"}]];
N27 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_28_2.mx"}]];
N28 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_29_2.mx"}]];
N29 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_30_2.mx"}]];
N30 = SelfEnergyFinite*kappa^2;


Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_3_2c.mx"}]];
N3c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_4_2c.mx"}]];
N4c = SelfEnergyFinite*kappa^2;


Nsum1b = makeFiniteAmplitude[N19+N20+N21+N22+N23+N24+N25+N26+N27+N28+N29+N30,-1,D];

Nsum1 = Nsum1b;

Nsum2b = makeFiniteAmplitude[N19+N20+N21+N22+N23+N24+N25+N26+N27+N28+N29+N30,-2,D];

Nsum2 = Nsum2b;

Nct1 = makeFiniteAmplitude[N3c+N4c,-1,D];
Nct2 = makeFiniteAmplitude[N3c+N4c,-2,D];

(* Particle 2 *)

SelfEnergyTotal = 0;



Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_26_2.mx"}]];
C26 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_27_2.mx"}]];
C27 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_28_2.mx"}]];
C28 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_29_2.mx"}]];
C29 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_30_2.mx"}]];
C30 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_31_2.mx"}]];
C31 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_32_2.mx"}]];
C32 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_33_2.mx"}]];
C33 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_34_2.mx"}]];
C34 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_36_2.mx"}]];
C36 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_37_2.mx"}]];
C37 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_38_2.mx"}]];
C38 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_39_2.mx"}]];
C39 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_40_2.mx"}]];
C40 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_41_2.mx"}]];
C41 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_42_2.mx"}]];
C42 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_43_2.mx"}]];
C43 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_44_2.mx"}]];
C44 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_45_2.mx"}]];
C45 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_46_2.mx"}]];
C46 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_47_2.mx"}]];
C47 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_48_2.mx"}]];
C48 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_49_2.mx"}]];
C49 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_50_2.mx"}]];
C50 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_51_2.mx"}]];
C51 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_52_2.mx"}]];
C52 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_53_2.mx"}]];
C53 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_56_2.mx"}]];
C56 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_57_2.mx"}]];
C57 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_58_2.mx"}]];
C58 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_59_2.mx"}]];
C59 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_60_2.mx"}]];
C60 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_61_2.mx"}]];
C61 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_62_2.mx"}]];
C62 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_63_2.mx"}]];
C63 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_64_2.mx"}]];
C64 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_65_2.mx"}]];
C65 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_66_2.mx"}]];
C66 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_67_2.mx"}]];
C67 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_69_2.mx"}]];
C69 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_70_2.mx"}]];
C70 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_71_2.mx"}]];
C71 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_72_2.mx"}]];
C72 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_73_2.mx"}]];
C73 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_74_2.mx"}]];
C74 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_4_2c.mx"}]];
C4c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_5_2c.mx"}]];
C5c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_6_2c.mx"}]];
C6c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_7_2c.mx"}]];
C7c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_8_2c.mx"}]];
C8c = SelfEnergyFinite*kappa^2;



Csum1b = makeFiniteAmplitude[C28+C29+C30+C31+C32+C33+C34+C36+C37+C38+C39+C40+C41+C42+C43+C44+C45+C46+C47+C48+C49+C50+C51+C52+C53+C56+C57+C58+C59+C60+C61+C62+C63+C64+C65+C66+C67+C69+C70+C71+C72+C73+C74,-1,D];
Csum1 = Csum1b;

Csum2b = makeFiniteAmplitude[C28+C29+C30+C31+C32+C33+C34+C36+C37+C38+C39+C40+C41+C42+C43+C44+C45+C46+C47+C48+C49+C50+C51+C52+C53+C56+C57+C58+C59+C60+C61+C62+C63+C64+C65+C66+C67+C69+C70+C71+C72+C73+C74,-2,D];
Csum2 = Csum2b;


Cct1 = makeFiniteAmplitude[C4c+C5c+C6c+C7c+C8c,-1,D];
Cct2 = makeFiniteAmplitude[C4c+C5c+C6c+C7c+C8c,-2,D];

{Csum1,Csum2,Nsum1,Nsum2,Cct1,Cct2,Nct1,Nct2}={Csum1,Csum2,Nsum1,Nsum2,Cct1,Cct2,Nct1,Nct2}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
/. MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]]/. MassBuilderP -> Momentum[p]/. MassBuilderEpsilon -> \[Epsilon]/. MassBuilderZeta -> \[Zeta]/. MassBuilderAe[0] -> 0\
/. D -> 4 - 2 \[Epsilon]/. Pair[p, p] -> p^2/. D -> 4 - 2 \[Epsilon]/. Momentum[p] -> p\
/. MassBuilderB[MChi, mz] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[mz, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[MChi, mw] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[mw, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[ma, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, ma}}]\
/. MassBuilderB[MChi, ma] -> TBI[4, MChi^2, {{1, MChi}, {1, ma}}]\
/. MassBuilderB[mf, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mf}}]\
/. MassBuilderB[MChi, mf] -> TBI[4, MChi^2, {{1, MChi}, {1, mf}}]\
/. MassBuilderB[mt, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mt}}]\
/. MassBuilderB[MChi, mt] -> TBI[4, MChi^2, {{1, MChi}, {1, mt}}]\
/.p->MChi/.STW->sw/.CTW->cw;
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];


Simplify[Csum2+Cct2 -Nsum2-Nct2/.mz->mw/cw]


Simplify[Nsum1+Nct1/.mz->mw/cw]


Simplify[Nct2/.mz->mw/cw]
Simplify[Nsum2/.mz->mw/cw]


Simplify[Cct2/.mz->mw/cw]
Simplify[Csum2/.mz->mw/cw]


CForm[FullSimplify[dZW1]]


CForm[FullSimplify[dMWsq1]]


CForm[FullSimplify[dMZsq1]]


CForm[FullSimplify[dZAA1]]


CForm[FullSimplify[dZZA1]]


(* ::InheritFromParent:: *)
(**)


CForm[FullSimplify[dZAZ1]]


CForm[FullSimplify[dZZZ1]]


FullSimplify[dZW1]


(* Check that counter-term couplings agree with Ibe et al. (2013) *)
g=e/sw;
gammaZ = -e^2(-5/3);
ZgammaZ = (-e*g/cw)*(11/6-5*sw^2/3);
ZZ= (-g^2/cw^2)( -11/6+11*sw^2/3-5 sw^4/3);
WZ=-g^2 * (-11/6);

ZgammaM = (-e*g/cw)*(2-2 sw^2)*mz^2;
ZM = (-g^2/cw^2)*(-1+6 sw^2- 4 sw^4)*mz^2;
WM = (-g^2)*(-1+2 sw^2)*mz^2;


FullSimplify[gammaZ- dZAA1]

FullSimplify[ZZ-dZZZ1/.mz->mw/cw]
FullSimplify[WZ-dZW1/.mz->mw/cw]
FullSimplify[ZgammaZ+(1/2)(dZZA1+dZAZ1)/.mz->mw/cw]

FullSimplify[dMZsq1+dZZZ1*mz^2 + ZM /.mz->mw/cw]

FullSimplify[dMWsq1+dZW1*mw^2 + WM /.mz->mw/cw]

FullSimplify[ZgammaM-(1/2)dZZA1*mz^2/.mz->mw/cw]




